function [X_,Y_,out] =  fcestimate(X_,Y_,iadj,params)
%FCESTIMATE Summary of this function goes here
%   Detailed explanation goes here
model = params.model;
optimopts = params.optimopts;
dispvec = X_-Y_;
x = dispvec(:,iadj);
% for non focus axis reject outliers based on vector norm. This should be
% (roughly) constant for non curvature directions
vcomp = dispvec(:,setdiff(1:3,iadj));
medvcomp = median(vcomp);
normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
normvcomp = sqrt(sum(normvcomp.*normvcomp,2));
if length(normvcomp)>20 & 0
    validinds = normvcomp<util.get1DThresh(normvcomp,20,.95);
else
    validinds = 1:length(normvcomp);
end
X_ = X_(validinds,:);
Y_ = Y_(validinds,:);
x = x(validinds,:);
%% rejection based on kdx
if iadj==1 % x-neighbor
    % x : x-displacement
    % y : y-location
    y = X_(:,2);
    bw = [2 220];
elseif iadj==2 % y-neighbor
    % x : y-displacement
    % y : x-location
    y = X_(:,1);
    bw=[3 100];
else % z-neighbor
    % x : z-displacement
    % y : y-location (not too much on x as cut is on y direction)
    y = X_(:,2);
    bw=[2 220];
end

if 1
    x_inline = x;
    y_inline = y;
else
    % build a probabilistic model of displacement vectors
    N = 101;
    gridx = linspace(min(x),max(x),N);
    gridy = linspace(min(y),max(y),N);
    [density,bw] = util.ksdensity2d([x y],gridx,gridy,bw);density=density'/max(density(:));
    [xmin,ix] = min(pdist2(x,gridx'),[],2);
    [ymin,iy] = min(pdist2(y,gridy'),[],2);
    idx = sub2ind([N,N],iy,ix);
    prob_inliers = density(idx)>max(density(idx))*.25;
    x_inline = x(prob_inliers,:);
    y_inline = y(prob_inliers,:);
    %%
    % fit curve model
    [~,im] = max(density,[],iadj);
    % arguably most important variable for successful fit
    sgn = 2*((max(im)==im(end) | max(im)==im(1))-.5);
end

if isfield(params,'init')
    pinit = params.init(iadj,:);
else
    pinit = [median(y) 1e-5 median(x)];
end
%%
% turn off warning
warning off
%TODO: very inefficient, find a robust way to pick inflection sign
[out1,r1] = nlinfit(y_inline, x_inline, model, pinit,optimopts);
[out2,r2] = nlinfit(y_inline, x_inline, model, pinit.*[1 -1 1],optimopts);
if norm(r1)<norm(r2);out = out1;else out = out2;end

if abs(out(2))<sqrt(eps) % mostlikely a line, p(1) is not reliable
    out(2) = 0; % to prevent scaling error in fc images
end
warning on
%%
% outlier rejection based on parametric model
xest = feval(model,out,y);
outliers = abs(x-xest)>2;
X_ = X_(~outliers,:);
Y_ = Y_(~outliers,:);

% %%
% figure, 
% plot(y_inline(~outliers),x_inline(~outliers),'.')
% hold on
% plot(y_inline,feval(model,out,y_inline),'.')
% 
% %%
% out=out1
% xest = feval(model,out,y);
% xline = feval(model,out,min(y):max(y));
% outliers = abs(x-xest)>2;
% % X_ = X_(~outliers,:);
% % Y_ = Y_(~outliers,:);
% 
% shift = median(dispvec);
% figure
% myplot3(X_,'.')
% hold on
% X_2 = X_;
% X_2(:,2) = xest;
% myplot3(X_2,'o')
% axis equal



end

