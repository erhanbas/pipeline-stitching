% applies descriptor/control pair parameters on a tile
% inputs:
%     tile ids: input tiles used in estimation
%     descriptor_file: descriptor match file, both in x/y/z
%     yml_file: transformation file used to map input tiles to um locations
addpath(genpath('./common'))
addpath(genpath('./visualization_functions'))
brain = '2018-08-01';
if 0
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
    % idx = 8569;
    ctrl = TileEstimator();
    ctrl.Vecfield = vecfield3D;
    ctrl.Scopeloc = scopeloc;
    ctrl.Neigs = neighbors;
    ctrl.pixres = pixres;
    ctrl.Regpts = regpts;
    ctrl.Paireddescriptor = paireddescriptor{1};
    ctrl.Scopeparams = scopeparams{1};
    
    %% mrse affine
    % stats: mean([xyz, x+1, x-1, y+1, y-1, z+1, z-1])
    numtile = length(vecfield3D.path);
    [Sest,Sres,Aest,Ares,Cres,residual_onx,residual_ony,residual_onz,stats] = deal(cell(1,numtile));
    tic
    parfor it = 1:numtile
        %%
        [Sest{it},Sres{it}] = ctrl.estimateStage(it);
        [Aest{it},Ares{it}] = ctrl.estimateAffine(it);
        [Cres{it},residual_onx{it}, residual_ony{it}, residual_onz{it}, stats{it}] = ctrl.estimateResidual4ctrl(it);
    end
    toc
    save ./visualization_functions/residual_init_results.mat -v7.3 ...
        Sest Sres Aest Ares Cres residual_onx residual_ony residual_onz stats
end
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% load ./visualization_functions/residual_init_results.mat Sest Sres Aest Ares Cres residual_onx residual_ony residual_onz stats
load ./visualization_functions/residual_iterated_results.mat Sest Sres Aest Ares Cres residual_onx residual_ony residual_onz stats


%%
% residual_onz_unbiased = residual_onz;
% for it=1:numtile
%     residual_onz_unbiased{it} = residual_onz_unbiased{it}-mean(residual_onz_unbiased{it});
% end

%%

%%
% load ./visualization_functions/residual_results.mat Sest Sres Aest Ares Cres residual_onx residual_ony residual_onz stats
%[S_mse,A_mse] = deal(nan(1,numtile));
[C_num, C_mse,C_msex,C_msey,C_msez] = deal(nan(1,numtile));
[CX_num, CX_mse,CX_msex,CX_msey,CX_msez] = deal(nan(1,numtile));
[CY_num, CY_mse,CY_msex,CY_msey,CY_msez] = deal(nan(1,numtile));
[CZ_num, CZ_mse,CZ_msex,CZ_msey,CZ_msez] = deal(nan(1,numtile));

mse_desc = nan(1,numtile);
for it = 1:numtile
    %inputval = inputval - mean(inputval);
    %if ~isnan(Aest{it}); S_mse(it) = ctrl.meanSqrt(Sres{it}); end
    %if ~isnan(Aest{it}); A_mse(it) = ctrl.meanSqrt(Ares{it}); end

    inputval = Cres{it};
    if ~isnan(Aest{it}) & ~isempty(residual_onz{it}) & ~isempty(inputval); [C_num(it), C_mse(it), C_msex(it), C_msey(it), C_msez(it)] = ctrl.tileStat(inputval); end
    inputval = residual_onx{it};
    if ~isnan(Aest{it}) & ~isempty(residual_onz{it}) & ~isempty(inputval); [CX_num(it), CX_mse(it), CX_msex(it), CX_msey(it), CX_msez(it)] = ctrl.tileStat(inputval); end
    inputval = residual_ony{it};
    if ~isnan(Aest{it}) & ~isempty(residual_onz{it}) & ~isempty(inputval); [CY_num(it), CY_mse(it), CY_msex(it), CY_msey(it), CY_msez(it)] = ctrl.tileStat(inputval); end
    inputval = residual_onz{it};
    if ~isnan(Aest{it}) & ~isempty(residual_onz{it}) & ~isempty(inputval); [CZ_num(it), CZ_mse(it), CZ_msex(it), CZ_msey(it), CZ_msez(it)] = ctrl.tileStat(inputval); end
end

% save mse_results_iterated C_num C_mse C_msex C_msey C_msez CX_num CX_mse CX_msex CX_msey CX_msez CY_num CY_mse CY_msex CY_msey CY_msez CZ_num CZ_mse CZ_msex CZ_msey CZ_msez scopeloc finterior -v7.3

%%
figure
h_ctrl = histogram(C_mse,'BinWidth',.025);

%%
% get mask/inds for interior. This is to prevent outliers due to gelatin
[interior] = valid_inds(scopeloc,1);
maskslices = scopeloc.gridix(:,3) < 1516 & scopeloc.gridix(:,3) > 1444;
finiteTiles = isfinite(C_mse(:));
finterior = find(interior & maskslices & finiteTiles);

%%
residual_stats = nan(numtile,8);
for it = 1:numtile
    if ~isnan(Aest{it}) & ~isempty(residual_onz{it})
        residual_stats(it,:) = [[size(Cres{it},1), size(residual_onx{it},1), size(residual_ony{it},1), size(residual_onz{it},1)], ...
            [ctrl.meanSqrt(Cres{it}), ctrl.meanSqrt(residual_onx{it}), ctrl.meanSqrt(residual_ony{it}-mean(residual_ony{it},'omitnan')), ctrl.meanSqrt(residual_onz{it}-mean(residual_onz{it},'omitnan'))]
            ];
    end
end
% save residual_stats residual_stats finterior
%%
figure
h_ctrl = histogram(residual_stats(finterior,5),'BinWidth',.025);



%%
figs = resStats();
figs.hist_mse(33,C_msex,C_msey,C_msez,finterior) %#ok<FNDSB>
title('residual')
% export_fig(fullfile('./visualization_figures','residual_iter_histogram_x_y_z.png'),'-transparent')

%%
finterior = find(interior & maskslices);
figs = resStats();
figs.hist_mse(34,CZ_msex,CZ_msey,CZ_msez,finterior) %#ok<FNDSB>
title('residual Z')
export_fig(fullfile('./visualization_figures','residualZ_iter_histogram_x_y_z.png'),'-transparent')

%%
save residual_results finterior C_mse C_msex C_msey C_msez CZ_mse CZ_msex CZ_msey CZ_msez

%%
[aa,bb]=max(C_mse(finterior));round([aa finterior(bb)])
it=finterior(bb);
round(scopeloc.loc(finterior(bb),:)*1e3)

%%



%%
figure(11), clf
% scatter3(round(scopeloc.gridix(:,1)),round(scopeloc.gridix(:,2)),round(scopeloc.gridix(:,3)),C_mse)
hold on
theseinds = 1:numtile;

scatter3(round(scopeloc.gridix(theseinds,1)),round(scopeloc.gridix(theseinds,2)),round(scopeloc.gridix(theseinds,3)),C_mse(theseinds))
zlabel('Z')
% axis equal

%%

[C_mse_sorted,sortinds] = sort(C_mse(finterior),'descend');
sortinds = finterior(sortinds);
disp(round(scopeloc.loc(sortinds(1:10),:)*1e3))
% disp(round(scopeloc.gridix(sortinds(1:10),1:3)))
testinds = sortinds(1:10);
sprintf('ind | C_mse | numDesc-x | numDesc-y | numDesc-z | Uniform | percent match')
round([testinds C_mse(testinds)' checkvalid(testinds,1:3) checkvalid(testinds,end-1)*100 checkvalid(testinds,end)*100])

% myplot3(round(scopeloc.gridix(:,1:3)),'o')
% myplot3(round(scopeloc.gridix(interior,1:3)),'r.')

% myplot3(round(scopeloc.gridix(interior,1:3)),'r.')


%%
clear resstats
resstats = resStats();
resstats.Aff = Aest;
resstats.Stg = Sest;
resstats.Stg_res = Sres;
resstats.Aff_res = Ares;
resstats.Ctrl_res = Cres;
resstats.Residual_onx = residual_onx;
resstats.Residual_ony = residual_ony;
resstats.Residual_onz = residual_onz;
resstats.init()
resstats.estimateMSE()

%%
%%
mean_res = nan(numtile,8);
for it = 1:numtile
    mean_res(it,:) = [[size(Cres{it},1), size(residual_onx{it},1), size(residual_ony{it},1), size(residual_onz{it},1)], ...
        [ctrl.meanSqrt(Cres{it}), ctrl.meanSqrt(residual_onx{it}), ctrl.meanSqrt(residual_ony{it}), ctrl.meanSqrt(residual_onz{it})]
        ];
end
% save mean_res_vs2 mean_res interior
%% viz
[aa,bb] = max(mean_res(finterior))
finterior(bb)
figure(11), cla
myplot3(round(scopeloc.gridix(:,1:3)),'o')
hold on
myplot3(round(scopeloc.gridix(interior,1:3)),'r.')
myplot3(round(scopeloc.gridix(finterior(bb),1:3)),'ks')
%% histograms
figs = resStats()
figs.hist_mse(32,C_mse,A_mse,S_mse,finterior)
% figs.boxplt(23,C_mse,A_mse,S_mse,finterior)
% export_fig(fullfile('./visualization_figures','residual_histogram.png'),'-transparent')
%% box plot

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










