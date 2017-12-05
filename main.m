function [outputArgs] = main()
%STICHING pipeline. Reads scope generated json file and returns a yml
%configuration file that goes into renderer. Requires calling cluster jobs
%to create subresults, i.e. descriptors. These functions can also run in
%local machines with proper settings.
%
% [OUTPUTARGS] = STICHING(jsonfile)
%
% Inputs:
%
% Outputs:
%
% Examples:
%
% Provide sample usage code here
%
% See also: List related files here

% $Author: base $	$Date: 2016/09/21 11:52:40 $	$Revision: 0.1 $
% Copyright: HHMI 2016

addpath(genpath('./common'))
addpath(genpath('./functions'))
brain = '2017-11-17';
tag='';

%% %%%%%%%%%%%%
directionMap = containers.Map({'-X','-Y','X','Y','-Z','Z'},[ 2, 3, 4, 5, 6, 7]);
directions = 'Z';

%% PATCH, fix this after dm11/tier2 merge
% raw input to descriptor generotion
inputfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/data/%s/Tiling',brain);
experimentfolder = sprintf('/nrs/mouselight/cluster/classifierOutputs/%s',brain);
% experimentfolder = sprintf('/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/%s%s/',brain,tag)

descriptorfolder = fullfile(experimentfolder,'classifier_output');
matfolder = fullfile(experimentfolder,'matfiles/');
scopefile = fullfile(matfolder,'scopeloc.mat');
desc_ch = {'1'};
descriptorfile = fullfile(matfolder,sprintf('descriptors_ch%s.mat',desc_ch{:})); % accumulated descriptor file
matchedfeatfile = fullfile(matfolder,sprintf('feats_ch%s.mat',desc_ch{:})); % accumulated descriptor file

mkdir(matfolder)
unix(sprintf('umask g+rxw %s',matfolder))
unix(sprintf('chmod g+rxw %s',matfolder))
%% 0: INTIALIZE
% read scope files and populate stage coordinates
if 0
    newdash = 1; % set this to 1 for datasets acquired after 160404
    [scopeloc] = getScopeCoordinates(inputfolder,newdash);% parse from acqusition files
    [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
    save(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
end

%% 1: LOAD MATCHED FEATS
if 0
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder');
    directions = 'Z';
    checkversion = 1; % 1: loads the version with "checkversion" extension and overwrites existing match if there are more matched points
    % load finished tile matches. find badly matched or missing tile pairs
    [regpts,featmap] = loadMatchedFeatures(scopeloc,descriptorfolder,directions,checkversion);

    save(fullfile(matfolder,'regpts.mat'),'regpts','featmap')
    if ~exist(fullfile(matfolder,'regpts_1stiter.mat'),'file')
        % faster to make a copy 
        unix(sprintf('cp %s %s',fullfile(matfolder,'regpts.mat'),fullfile(matfolder,'regpts_1stiter.mat')))
        % save(fullfile(matfolder,'regpts_1stiter.mat'),'regpts','featmap')
    end
end

%%
if 1 % iterate on missing tiles
    addpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'),'-end')
    runlocal = 1;
    pointmatch_task(brain,runlocal)
    rmpath(genpath('/groups/mousebrainmicro/home/base/CODE/MATLAB/pipeline/zmatch_pipe'))
end

%% 2 scope params estimation
%%
if 1
    %%
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    scopeacqparams = util.readScopeFile(fileparts(scopeloc.filepath{1}));
    params.scopeacqparams = scopeacqparams;
    params.viz = 0;
    params.debug = 0;
    params.Ndivs = 4;
    params.Nlayer = 4;
    params.htop = 5;
    params.expensionratio = 1;
    params.imagesize = [1024 1536 251];
    params.imsize_um = [scopeacqparams.x_size_um scopeacqparams.y_size_um scopeacqparams.z_size_um];
    params.overlap_um = [scopeacqparams.x_overlap_um scopeacqparams.y_overlap_um scopeacqparams.z_overlap_um];
    params.order = 1;
    params.applyFC = 1;
    params.beadparams = [];%beadparams;

    descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);
    [paireddescriptor, scopeparams, R, curvemodel,scopeparams_, paireddescriptor_,curvemodel_] = ...
        estimateFCpertile(descriptors,neighbors,scopeloc,params); %#ok<ASGLU>
    %%
    save(descriptorfile,'descriptors','-v7.3')
    save(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
        'scopeparams', 'R', 'curvemodel','scopeparams_', 'paireddescriptor_', ...
        'curvemodel_','params')
end

%%
if 1
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'scopeparams_pertile'),'scopeparams')
    load(fullfile(matfolder,'regpts'),'regpts')
    videofile = sprintf('./videos/%s-1stiter-ch1',brain)
    % descriptorMatchQuality(regpts,scopeparams,scopeloc,videofile)
%     createThumb(regpts,scopeparams,scopeloc,videofile)
    descriptorMatchQualityHeatMap(regpts,scopeparams,scopeloc,videofile)
end

%%
if 1
    load(scopefile,'scopeloc','neighbors','experimentfolder','inputfolder')
    load(fullfile(matfolder,'regpts'),'regpts')
    load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
        'scopeparams', 'R', 'curvemodel','scopeparams_', 'paireddescriptor_', ...
        'curvemodel_','params')
    
    vecfield3D = vectorField3D(params,scopeloc,regpts,scopeparams_,curvemodel_,[]);
    if 1
        save(fullfile(matfolder,sprintf('%s_%s',datestr(now,'mmddyyHHMMSS'),'vecfield3D')),'vecfield3D','params')
        save(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
    end
end

% 4
load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
load(fullfile(matfolder,'vecfield3D'),'vecfield3D','params')
vecfield = vecfield3D;

% checkthese = [1 4 5 7]; % 0 - right - bottom - below
% [neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
params.big = 1;
params.ymldims = [params.imagesize 2];%[1024 1536 251 2]
sub = 0;
params.root = vecfield.root;

if sub
    targetidx = getTargetIDx(scopeloc,neighbors);
    params.outfile = fullfile(experimentfolder,sprintf('%s_sub.control.yml',date));
else
    targetidx = 1:size(scopeloc.gridix,1);
    params.outfile = fullfile(experimentfolder,sprintf('%s.control.yml',date));
end

writeYML(params, targetidx(:)', vecfield)
unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache.yml')))
%
if ~sub
    params.big=0
    params.outfile = sprintf('%s/%s.old.control.yml',experimentfolder,date);
    writeYML(params, targetidx(:)', vecfield)
    unix(sprintf('cp %s %s',params.outfile,fullfile(experimentfolder,'tilebase.cache_old.yml')))
end
return
%% 0.1: FLAT RUN
% generate yml for without any optimization
if 0
    %%
    load(scopefile,'scopeloc','neighbors','imsize_um','experimentfolder','inputfolder')
    if scope==1
        scope1_beadparams = load('./beadparams/scope1_beadparams');
        scopeparams = scope1_beadparams.scope1_beadparams;
    else
        scope2_beadparams = load('./beadparams/scope2_beadparams');
        scopeparams = scope2_beadparams.scope2_beadparams;
    end
    vecfield = vectorField_flatrun(params,scopeloc,scopeparams,2);

    load ./matfiles/xypaireddescriptor paireddescriptor R curvemodel
    [scopeparams,scopeparams_,paireddescriptor_,curvemodel_] = homographyPerTile6Neighbor(...
        beadparams,neighbors,scopeloc,paireddescriptor,R,curvemodel,imsize_um);
    vecfield3D_flat_4neig = vectorField_flatrun_pertile(params,scopeloc,scopeparams_,curvemodel_,[]);
    save pertile_4neig scopeparams scopeparams_ paireddescriptor_ curvemodel_ vecfield3D_flat_4neig
end
%% stitching quality test
if 0
    load(fullfile(matfolder,'scopeloc'),'scopeloc','imsize_um','experimentfolder','inputfolder')
    load(fullfile(matfolder,'vecfield'),'vecfield','params')
    %%
    clc
    params.big = 1;
    params.dims = [params.imagesize 2]%[1024 1536 251 2]
    sub = 0;
    inds_ = inds(1)';
    neigs = neighbors(inds_,checkthese);
    targetidx = neigs([1 3])
    params.root = vecfield.root;
    %%
    %
    params.outfile = sprintf('%s%s-%d_%d_sub_1.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(1)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y1-ccx2-fixxed')
    %%
    
    params.outfile = sprintf('%s%s-%d_%d_sub_2.tilebase.cache.yml',experimentfolder,date,targetidx);
    params.outfile
    vecfield_ = vecfield;
    vecfield_.path{targetidx(2)} = '/00000';
    writeYML(params, targetidx(:)', vecfield_)
    paramoutput = '/groups/mousebrainmicro/mousebrainmicro/cluster/Stitching/2017-02-13/set_parameters_sub'
    converter(params.outfile,paramoutput,'y2-ccx2-fixxed')
end
end
