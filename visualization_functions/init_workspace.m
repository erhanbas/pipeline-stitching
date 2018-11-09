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

save ../visualization_functions/residual_init_results.mat -v7.3 ...
    Sest Sres Aest Ares Cres residual_onx residual_ony residual_onz stats

disp('DONE')