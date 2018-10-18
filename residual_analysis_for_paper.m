% applies descriptor/control pair parameters on a tile
% inputs: 
%     tile ids: input tiles used in estimation
%     descriptor_file: descriptor match file, both in x/y/z
%     yml_file: transformation file used to map input tiles to um locations

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

%% estimate the residule
idx_center = 5000;
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
% baseline: apply stage transformation
res_1 = ctrl.residualStageOnZ(idx_center);
res_2 = ctrl.residualAffineFCOnZ(idx_center);


%%
tile_loc_um_center = scopeloc.loc(idx_center,:)*1e3; %in um
tile_loc_um_target = scopeloc.loc(idx_target,:)*1e3; %in um

tile_descs_center_um = tile_descs_center.*pixres + tile_loc_um_center;
tile_descs_target_um = tile_descs_target.*pixres + tile_loc_um_target;

% tile_loc_um_center - tile_loc_um_target
residual_1 = tile_descs_center_um - tile_descs_target_um;

%%
aff_center = ctrl.Vecfield.afftile(:,:,idx_center)/1e3; %in um
aff_target = ctrl.Vecfield.afftile(:,:,idx_target)/1e3; %in um

tile_descs_center_aff_um = aff_center * [tile_descs_center ones(size(tile_descs_center,1),1)]';
tile_descs_target_aff_um = aff_target * [tile_descs_target ones(size(tile_descs_target,1),1)]';

% tile_loc_um_center - tile_loc_um_target
residual_2 = [tile_descs_center_aff_um - tile_descs_target_aff_um]';

%% control points

aff_center = ctrl.Vecfield.afftile(:,:,idx_center)/1e3; %in um
aff_target = ctrl.Vecfield.afftile(:,:,idx_target)/1e3; %in um

tile_descs_center_aff_um = aff_center * [tile_descs_center ones(size(tile_descs_center,1),1)]';
tile_descs_target_aff_um = aff_target * [tile_descs_target ones(size(tile_descs_target,1),1)]';

% tile_loc_um_center - tile_loc_um_target
residual_2 = [tile_descs_center_aff_um - tile_descs_target_aff_um]';

%%
%estimated transform
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
figure(12), cla
myplot3(tile_descs_center,'+')
hold on
myplot3(xyz_center,'o')

% myplot3(bottom_tile_descs,'r+')










