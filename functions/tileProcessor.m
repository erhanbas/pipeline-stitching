function [descriptors,paireddescriptor_,curvemodel_,scopeparams_,paireddescriptor,curvemodel,scopeparams] = tileProcessor(scopeloc,descriptorfolder,desc_ch)

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

%% 
descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);

%%
[paireddescriptor,R,curvemodel] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params);

%%
[paireddescriptor_,curvemodel_] = match.curvatureOutlierElimination(paireddescriptor,curvemodel,scopeloc);
%%
% joint affine estimation
[scopeparams] = match.estimatejointaffine(paireddescriptor_,neighbors,scopeloc,params,curvemodel_,0);
%%
scopeparams_ = match.affineOutlierElimination( scopeparams );



