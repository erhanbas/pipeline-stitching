function [X_,Y_] = descriptorMatch4XY(X,Y,params)
% descriptor match for XY directions. Uses tighter constraints as there is
% no physical cut and initialization can be reliably estimated from beads

opt = params.opt;
projectionThr = params.projectionThr;
%%
% eliminate isolated descriptors
pDori = pdist2(X,Y);
[YnearestD,bb1]=min(pDori,[],1);
[XnearestD,bb2]=min(pDori,[],2);
% for X/Y
X = X(XnearestD<projectionThr,:);
Y = Y(YnearestD<projectionThr,:);
if size(X,1)<3 | size(Y,1)<3; [X_,Y_] = deal([]);return;end
    
%%
% initial match based on point drift
[Transform, C]=cpd_register(X,Y,opt);
%% check if match is found
pD = pdist2(X,Transform.Y);
[aa1,bb1]=min(pD,[],1);
[aa2,bb2]=min(pD,[],2);
keeptheseY = find([1:length(bb1)]'==bb2(bb1));
keeptheseX = bb1(keeptheseY)';

disttrim = aa1(keeptheseY)'<projectionThr;
X_ = X(keeptheseX(disttrim),:);
Y_ = Y(keeptheseY(disttrim),:);

