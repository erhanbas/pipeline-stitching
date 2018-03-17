function [descriptors,paireddescriptor,curvemodel,scopeparams] = tileProcessor_debug(scopeloc,descriptorfolder,desc_ch,params)

checkthese = [1 4 5 7]; % 0 - right - bottom - below
[neighbors] = buildNeighbor(scopeloc.gridix(:,1:3)); %[id -x -y +x +y -z +z] format
% accumulator for steps
[paireddescriptor,R,curvemodel,scopeparams]=deal([]);

%% get tile descriptors
if 0
    descriptors = getDescriptorsPerFolder(descriptorfolder,scopeloc,desc_ch);
else
    load('/nrs/mouselight/cluster/classifierOutputs/2017-09-25/matfiles/descriptors_ch1_.mat','descriptors')
end

%% descriptor match
% [paireddescriptor{end+1},R{end+1},curvemodel{end+1}] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params);
[paireddescriptor{1},R{1},curvemodel{1}] = match.xymatch(descriptors,neighbors(:,checkthese),scopeloc,params);

%%
% [paireddescriptor{end+1},curvemodel{end+1},unreliable] = match.curvatureOutlierElimination(paireddescriptor{end},curvemodel{end},scopeloc);
[paireddescriptor{2},curvemodel{2},unreliable] = match.curvatureOutlierElimination(paireddescriptor{1},curvemodel{1},scopeloc);
%%
% tile base affine
if params.singleTile
%     [scopeparams{1},scopeparams{2},paireddescriptor{end+1},curvemodel{end+1}] = homographyPerTile6Neighbor(...
%         params,neighbors,scopeloc,paireddescriptor{end},R,curvemodel{end});
    % debug
    [scopeparams{1},scopeparams{2},paireddescriptor{3},curvemodel{3}] = homographyPerTile6Neighbor(...
        params,neighbors,scopeloc,paireddescriptor{2},R,curvemodel{2});
    
else
    % joint affine estimation
    [scopeparams{end+1}] = match.estimatejointaffine(paireddescriptor{end},neighbors,scopeloc,params,curvemodel{end},0);
    [scopeparams{end+1}, paireddescriptor{end+1}, curvemodel{end+1}] = match.affineOutlierElimination( scopeloc, scopeparams{end}, paireddescriptor{end},curvemodel{end},unreliable );
end

% %%
%         matfolder='/nrs/mouselight/cluster/classifierOutputs/2017-09-25/matfiles/'
%         load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
%             'scopeparams', 'curvemodel','params')
%         %%
%         save(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
%             'scopeparams', 'curvemodel','params','-v7.3')
% 


