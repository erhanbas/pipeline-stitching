% applies descriptor/control pair parameters on a tile
% inputs: 
%     tile ids: input tiles used in estimation
%     descriptor_file: descriptor match file, both in x/y/z
%     yml_file: transformation file used to map input tiles to um locations

addpath(genpath('./common'))
addpath(genpath('./visualization_functions'))
brain = '2018-08-01';
desc_ch = {'0'};

experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file

load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
load(fullfile(matfolder,'scopeparams_pertile'),'scopeparams')
load(fullfile(matfolder,'regpts'),'regpts')
load(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
load(descriptorfile,'descriptors')
load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
    'scopeparams', 'curvemodel','params')

checkthese = [1 4 5 7]; % 0 - right - bottom - below
imsize_um = params.imsize_um;
FOV = [params.scopeacqparams.fov_x_size_um params.scopeacqparams.fov_y_size_um params.scopeacqparams.fov_z_size_um];
pixres = FOV./(params.imagesize-1);

input_tile_id = 1000;
neighboring_tile_id = neighbors(input_tile_id,:);
if 0
    [paireddescriptor,R,curvemodel] = xymatch_old(...
        descriptors,neighbors(:,checkthese),scopeloc,params);
    [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
        params,neighbors,scopeloc,paireddescriptor,R,curvemodel);
end
%% compile descriptor matches

% tot_num_tiles = length(scopeloc.filepath);
% pairGraph = cell(tot_num_tiles,tot_num_tiles);
% %%
% descset = paireddescriptor{3};
% for idx_tile = 1: tot_num_tiles
%     neigs = neighbors(idx_tile,:); %[id -x -y +x +y -z +z]
%     % onx
%     ix = idx_tile;
%     X = descset{ix}.onx.X;
%     Y = descset{ix}.onx.Y;
%     
%     iy = neigs(2);
%     X = descset{iy}.onx.X;
%     Y = descset{iy}.onx.Y;
% 
%     iy = neigs(5);
%     pairGraph{ix,iy}.ony = paireddescriptor{1}{ix}.ony;
%     regpts{ix}
%     pairGraph{ix,iy}.onz = paireddescriptor{1}{ix}.onz;
% 
%     descpairs
%     
%     
%     
% end
%% estimate the residule
idx = 8569;
ctrl = TileEstimator();
ctrl.Vecfield = vecfield3D;
ctrl.Scopeloc = scopeloc;
ctrl.Neigs = neighbors;
ctrl.pixres = pixres;
ctrl.Regpts = regpts;
ctrl.Paireddescriptor = paireddescriptor{end};
ctrl.Scopeparams = scopeparams{end};

%% get mask/inds for interior. This is to prevent outliers due to gelatin 
[interior] = valid_inds(scopeloc,1);

%% mrse affine
% stats: mean([xyz, x+1, x-1, y+1, y-1, z+1, z-1])
numtile = length(vecfield3D.path);
[Sest,Sres,Aest,Ares,Cres,residual_onx,residual_ony,residual_onz,stats] = deal(cell(1,numtile));
tic
parfor it = 1:numtile
    [Sest{it},Sres{it}] = ctrl.estimateStage(it);
    [Aest{it},Ares{it}] = ctrl.estimateAffine(it);
    [Cres{it},residual_onx{it}, residual_ony{it}, residual_onz{it}, stats{it}] = ctrl.estimateResidual4ctrl(it);
end
toc
% save workspace.mat -v7.3

%%
[S_mse,A_mse,C_mse] = deal(nan(1,numtile));
mse_desc = nan(1,numtile);
for it = 1:numtile
    if ~isnan(Aest{it}); S_mse(it) = mean(sqrt(sum(Sres{it}.^2,2))); end
    if ~isnan(Aest{it}); A_mse(it) = mean(sqrt(sum(Ares{it}.^2,2))); end
    if ~isnan(Aest{it}); C_mse(it) = mean(sqrt(sum(Cres{it}.^2,2)),'omitnan'); end
end

mean_res = nan(numtile,4);
for it = 1:numtile
    mean_res(it,:) = [ctrl.meanSqrt(Cres{it}),ctrl.meanSqrt(residual_onx{it}) ctrl.meanSqrt(residual_ony{it}) ctrl.meanSqrt(residual_onz{it})];
end

%% viz
finterior = find(interior);
[aa,bb] = max(mean_res(finterior))
finterior(bb)
figure(11), cla
myplot3(round(scopeloc.gridix(:,1:3)),'o')
hold on
myplot3(round(scopeloc.gridix(interior,1:3)),'r.')
myplot3(round(scopeloc.gridix(finterior(bb),1:3)),'ks')
%% histograms
% close all
clc
figure(123), clf

ax = gca;
h_ctrl = histogram(C_mse(finterior));
hold on
h_aff = histogram(A_mse(finterior));
h_st = histogram(S_mse(finterior));


C_st = [mean(C_mse(finterior),'omitnan') std(C_mse(finterior),'omitnan')];
A_st = [mean(A_mse(finterior),'omitnan') std(A_mse(finterior),'omitnan')];
S_st = [mean(S_mse(finterior),'omitnan') std(S_mse(finterior),'omitnan')];
leg1 = sprintf('Control points: %1.3f \\pm %1.3f',C_st(1),C_st(2));
leg2 = sprintf('Affine: %1.3f \\pm %1.3f',A_st(1),A_st(2));
leg3 = sprintf('Stage: %1.3f \\pm %1.3f',S_st(1),S_st(2));

leg = legend(leg1,leg2,leg3);
% legend('test \pm')
leg.FontSize = 32;
h_ctrl.BinWidth = .2;
h_aff.BinWidth = .2;
h_st.BinWidth = .2;
xlim([0 20])

title('Residual out of 11657 tiles', 'FontSize', 30)
xlabel('Residual magnitude in \mum', 'FontSize', 30)
ax.YAxis.FontSize = 24;
ax.XAxis.FontSize = 24;
export_fig(fullfile('./visualization_figures','residual_histogram.png'),'-transparent')
%% box plot
len = length(finterior);
C_fint = C_mse(finterior);
S_fint = S_mse(finterior);
A_fint = A_mse(finterior);
X_fint = [C_fint,A_fint,S_fint];

Origin = [repmat({'Ctrl'},len,1);repmat({'Affine'},len,1);repmat({'Stage'},len,1)];
leg1 = sprintf('%1.3f \\pm %1.3f',C_st(1),C_st(2));
leg2 = sprintf('%1.3f \\pm %1.3f',A_st(1),A_st(2));
leg3 = sprintf('%1.3f \\pm %1.3f',S_st(1),S_st(2));

figure(34), 
clf,
ax = cla;
hold on
boxplot(X_fint,Origin)
text(1-.1,median(C_fint,'omitnan'), leg1,'FontSize',12,'FontWeight','bold')
text(2-.1,median(A_fint,'omitnan'), leg2,'FontSize',12,'FontWeight','bold')
text(3-.1,median(S_fint,'omitnan'), leg3,'FontSize',12,'FontWeight','bold')

ylim([-1 30])
ylabel('Residual in \mum')

ax.YAxis.FontSize = 34;
ax.XAxis.FontSize = 44;
set(ax,'Ytick',0:5:32)
export_fig(fullfile('./visualization_figures','residual_boxplot.png'),'-transparent')

%%
% figure, 
% hist(mean_res(finterior),100)
% figure, hist(mean_res(finterior),100)
%%
disp(sprintf('MeanSE of residuals for %d tiles: %s %s %s',...
    sum(isfinite(A_mse)),...
    'based on stage |', 'based on affine |', 'based on ctrl'))

disp(sprintf('MeanSE of residuals for %d: %f | %f | %f',...
    sum(isfinite(A_mse)),...
    mean(S_mse,'omitnan'), mean(A_mse(isfinite(A_mse)),'omitnan'), mean(C_mse(isfinite(A_mse)),'omitnan')))

disp(sprintf('MedianSE of residuals for %d: %f | %f | %f',...
    sum(isfinite(A_mse)),...
    median(S_mse,'omitnan'), median(A_mse,'omitnan'), median(C_mse,'omitnan')))

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

tt = scopeloc.gridix(these_inds,1:3);
% close all
%%
% get indicies of 

figure, 
imshow(sum(tileImage,3),[])
figure, 
imshow(sum(out,3),[])
%%
[A_mse_sorted,A_mse_inds] = sort(A_mse(isfinite(A_mse)),'descend');
[C_mse_sorted,C_mse_inds] = sort(C_mse(isfinite(C_mse)),'descend');

bb = find(A_mse<C_mse);
aa = find(A_mse>C_mse);
cc = find(C_mse>30);

%
figure(13), cla
myplot3(round(scopeloc.loc(bb,:)*1e3),'o')
hold on
myplot3(round(scopeloc.loc(aa,:)*1e3),'r.')
myplot3(round(scopeloc.loc(cc,:)*1e3),'gs')

%%
ctrl.estimateResidual4ctrl(idx_center)



%%
[accumulator_stage,accumulator_descriptor] = ctrl.getDifferences(idx);

%%
idx_center = 8569;
[Aest_idx,res] = ctrl.estimateAffine(idx_center);
mean(sqrt(sum(res.^2,2)))

idx_target = neighbors(idx_center,end);
[Aest_idy,res_idy] = ctrl.estimateAffine(idx_target);

%
tile_descs_center = regpts{idx_center}.X;
tile_descs_target = regpts{idx_center}.Y;

x1 = tile_descs_center * Aest_idx + ctrl.Scopeloc.loc(idx_center,:)*1e3;
x2 = tile_descs_target * Aest_idy + ctrl.Scopeloc.loc(idx_target,:)*1e3;

x1-x2

%%
figure, 
myplot3(tile_descs_center,'o')
hold on
myplot3(tile_descs_target,'*')

%%
figure, 
myplot3(x1,'o')
hold on
myplot3(x2,'*')

%% estimate the residule
idx_center = 8569;
idx_target = neighbors(idx_center,end);

% along z direction
tile_descs_center = regpts{idx_center}.X;
tile_descs_target = regpts{idx_center}.Y;

%
ctrl = ControlPoints();
ctrl.Vecfield = vecfield3D;
ctrl.Scopeloc = scopeloc;
ctrl.Neigs = neighbors;
ctrl.pixres = pixres;
ctrl.Regpts = regpts;

xyz_center = ctrl.gridLocations(idx_center);
xyz_target = ctrl.gridLocations(idx_target);
control_points_center = ctrl.Vecfield.control(:,:,idx_center)/1e3; % in um
control_points_target = ctrl.Vecfield.control(:,:,idx_target)/1e3; % in um
%%
edgelist = drawPlates(xyz_center,[],[5 5 4]);

figure(111), cla
bbox = round([min([control_points_center;control_points_target]) max([control_points_center;control_points_target])]);
drawPlates(control_points_center,edgelist,[0 0 1],[5 5 4]);
drawPlates(control_points_target,edgelist,[1 1 0],[5 5 4]);
xlim([bbox(1) bbox(4)] + [-1 1]*10)
ylim([bbox(2) bbox(5)] + [-1 1]*10)
zlim([bbox(3) bbox(6)] + [-1 1]*10)
view([-27 15]),
axis off
% export_fig(fullfile('./visualization_figures','ctrl_top_tile.png'),'-transparent')
%
% add descriptor points with transforms
tile_descs_center_ctrl_um = ctrl.interpolatedInstance(xyz_center,control_points_center,tile_descs_center);
tile_descs_target_ctrl_um = ctrl.interpolatedInstance(xyz_target,control_points_target,tile_descs_target);
clc
myplot3(tile_descs_center_ctrl_um,{'bo','MarkerFaceColor','b'})
myplot3(tile_descs_target_ctrl_um,{'yo','MarkerFaceColor','y'})
%%
valid_center = any(~isnan(tile_descs_center_ctrl_um),2);
valid_target = any(~isnan(tile_descs_target_ctrl_um),2);
valid_descs = valid_center & valid_target;

% fix any stage shift problems
cent_stg_um_corr = cent_stg_um+(mean(target_stg_um)-mean(cent_stg_um));
cent_aff_um_corr = cent_aff_um+(mean(target_aff_um)-mean(cent_aff_um)); % need to add x&y
cent_ctrl_um_corr = cent_ctrl_um+(mean(target_ctrl_um,'omitnan')-mean(cent_ctrl_um,'omitnan'));



dif_stg_um = cent_stg_um_corr-target_stg_um;
dif_aff_um = cent_aff_um_corr-target_aff_um;
dif_ctrl_um = cent_ctrl_um_corr-target_ctrl_um;
max(abs(dif_stg_um))
max(abs(dif_aff_um))

%%
% baseline: apply stage transformation
[res_1, cent_stg_um, target_stg_um] = ctrl.residualStageOnZ(idx_center);
[res_2, cent_aff_um, target_aff_um] = ctrl.residualAffineFCOnZ(idx_center);
[res_3, cent_ctrl_um, target_ctrl_um] = ctrl.residualCtrlOnZ(idx_center);

% if we were to rely on just the stage, residual will be
err_stage_um = (norm(res_1)); % in um
% correct the shift/translation only based on match
err_stage_translation_um = sqrt(sum((res_1-mean(res_1)).^2,2)); % in um
% if we were to rely on just the stage, residual will be
err_aff_um = (norm(res_2)); % in um
% correct the shift/translation only based on match
err_aff_translation_um = sqrt(sum((res_2-mean(res_2)).^2,2)); % in um
%%
mean(err_stage_translation_um)
mean(err_aff_translation_um)

% mean(sqrt(sum(res_3.^2,2)),'omitnan')
%%

err_aff_um_1 = norm(res_2(:,1:3));
err_aff_um_2 = norm(res_2(:,4:6));
% correct the affine then shift/translation only based on match
err_aff_translation_um = (norm(res_2(:,4:6)-mean(res_2(:,4:6)))); % in um
%%
edgelist = drawPlates(xyz_center,[],[5 5 4]);
close all

figure(111), cla
drawPlates(control_points_center,edgelist,[0 0 1],[5 5 4]);
drawPlates(control_points_target,edgelist,[1 1 0],[5 5 4]);

%%
figure(11), cla,
set(gca,'Zdir','reverse')
myplot3(cent_stg_um+(mean(target_stg_um)-mean(cent_stg_um)),'+')
hold on
myplot3(target_stg_um,'o')
axis tight equal
title('STAGE only')

figure(12), cla,
set(gca,'Zdir','reverse')
myplot3(cent_aff_um+(mean(target_aff_um)-mean(cent_aff_um)),'+')
hold on
myplot3(target_aff_um,'o')
axis tight equal
title('Affine only')

figure(13), cla,
set(gca,'Zdir','reverse')
myplot3(cent_ctrl_um+(mean(target_ctrl_um)-mean(cent_ctrl_um)),'+')
hold on
myplot3(target_ctrl_um,'o')
axis tight equal
title('ctrl')

%%
%estimated transform
figure(121), cla
myplot3(tile_descs_center-mean(tile_descs_center),'bo')
hold on
myplot3(tile_descs_target-mean(tile_descs_target),'r+')
axis tight equal
% myplot3(xyz_center,'rs')
%%
figure(13), cla,
myplot3(xyz_center,'+')
hold on
myplot3(tile_descs_center,'o')
set(gca,'Zdir','reverse')

% %%
% tile_loc_um_center = scopeloc.loc(idx_center,:)*1e3; %in um
% tile_loc_um_target = scopeloc.loc(idx_target,:)*1e3; %in um
% 
% tile_descs_center_um = tile_descs_center.*pixres + tile_loc_um_center;
% tile_descs_target_um = tile_descs_target.*pixres + tile_loc_um_target;
% 
% % tile_loc_um_center - tile_loc_um_target
% residual_1 = tile_descs_center_um - tile_descs_target_um;
% 
% %%
% aff_center = ctrl.Vecfield.afftile(:,:,idx_center)/1e3; %in um
% aff_target = ctrl.Vecfield.afftile(:,:,idx_target)/1e3; %in um
% 
% tile_descs_center_aff_um = aff_center * [tile_descs_center ones(size(tile_descs_center,1),1)]';
% tile_descs_target_aff_um = aff_target * [tile_descs_target ones(size(tile_descs_target,1),1)]';
% 
% % tile_loc_um_center - tile_loc_um_target
% residual_2 = [tile_descs_center_aff_um - tile_descs_target_aff_um]';

%% control points
FxU = scatteredInterpolant(xyz_center,control_points_center(:,1),'linear','nearest');
FxV = scatteredInterpolant(xyz_center,control_points_center(:,2),'linear','nearest');
FxW = scatteredInterpolant(xyz_center,control_points_center(:,3),'linear','nearest');
tile_descs_center_ctrl_um = [FxU(tile_descs_center) FxV(tile_descs_center) FxW(tile_descs_center)];

FxU = scatteredInterpolant(xyz_target,control_points_target(:,1),'linear','nearest');
FxV = scatteredInterpolant(xyz_target,control_points_target(:,2),'linear','nearest');
FxW = scatteredInterpolant(xyz_target,control_points_target(:,3),'linear','nearest');
tile_descs_target_ctrl_um = [FxU(tile_descs_target) FxV(tile_descs_target) FxW(tile_descs_target)];

residual_3 = tile_descs_center_ctrl_um - tile_descs_target_ctrl_um;
[res_3, cent_um, target_um] = ctrl.residualCtrlOnZ(idx_center);
%%
%estimated transform
figure(15), cla,
myplot3(control_points_center,'bo')
hold on
myplot3(tile_descs_center_ctrl_um,'m+')

% myplot3(control_points_target,'go')
myplot3(tile_descs_target_ctrl_um,'rd')

set(gca,'Zdir','reverse')
%%
figure(14), cla,
myplot3(xyz_center,'+')
set(gca,'Zdir','reverse')


%%
figure(14), cla
myplot3(tile_descs_center_aff_um','+')
hold on
myplot3(tile_descs_target_aff_um','o')

%% for x&y 
% apply affine without FC correction


% apply affine with FC correction


% apply affine with FC/homography estimation


%%
%estimated transform

% myplot3(bottom_tile_descs,'r+')










