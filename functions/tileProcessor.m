function [descriptors,paireddescriptor,curvemodel,scopeparams] = tileProcessor(scopeloc,descriptorfolder,desc_ch)

% paramater setting for descrtiptor match
scopeacqparams = util.readScopeFile(fileparts(scopeloc.filepath{1}));
% params.sample = brain;
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
params.beadparams = [];%PLACEHOLDER FOR BEADS, very unlikely to have it...

checkthese = [1 4 5 7]; % 0 - right - bottom - below
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% accumulator for steps
[paireddescriptor,R,curvemodel,scopeparams]=deal([]);

%% get tile descriptors
descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);

%% descriptor match
[paireddescriptor{end+1},R{end+1},curvemodel{end+1}] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params);
%%
[paireddescriptor{end+1},curvemodel{end+1},unreliable] = match.curvatureOutlierElimination(paireddescriptor{end},curvemodel{end},scopeloc);
%%
% joint affine estimation
[scopeparams{end+1}] = match.estimatejointaffine(paireddescriptor{end},neighbors,scopeloc,params,curvemodel{end},0);
%%
[scopeparams{end+1}, paireddescriptor{end+1}, curvemodel{end+1}] = match.affineOutlierElimination( scopeloc, scopeparams{end}, paireddescriptor{end},curvemodel{end},unreliable );
%%
% %%
%         matfolder='/nrs/mouselight/cluster/classifierOutputs/2017-09-25/matfiles/'
%         load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
%             'scopeparams', 'R', 'curvemodel', 'paireddescriptor_', ...
%             'curvemodel_','params')
%         %%
%         save(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
%             'scopeparams', 'R', 'curvemodel', 'paireddescriptor_', ...
%             'curvemodel_','params','-v7.3')



