function [X_,Y_,out,valid] =  fcestimate(X_,Y_,iadj,params)
%FCESTIMATE Summary of this function goes here
%   Detailed explanation goes here
valid = 0;
viz = params.viz;
% model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
model = params.model;
optimopts = params.optimopts;
dispvec = X_-Y_;
y = dispvec(:,iadj);

% for non focus axis reject outliers based on vector norm. This should be
% (roughly) constant for non curvature directions
if 1
    validinds = 1:length(y);
elseif length(normvcomp)>20 & 0
    vcomp = dispvec(:,setdiff(1:3,iadj));
    medvcomp = median(vcomp);
    normvcomp = vcomp-ones(size(vcomp,1),1)*medvcomp;
    normvcomp = sqrt(sum(normvcomp.*normvcomp,2));
    validinds = normvcomp<util.get1DThresh(normvcomp,20,.95);
else 0
end

X_ = X_(validinds,:);
Y_ = Y_(validinds,:);
y = y(validinds,:);

x = X_(:,setdiff([1 2],iadj));

if isfield(params,'init')
    pinit = params.init(iadj,:);
else
    pinit = [median(y) 1e-5 median(x)];
end

%%
% turn off warning
warning off
%TODO: very inefficient, find a robust way to pick inflection sign
[out1,r1] = nlinfit(x, y, model, pinit,optimopts);
[out2,r2] = nlinfit(x, y, model, pinit.*[1 -1 1],optimopts);
if norm(r1)<norm(r2);out = out1;else out = out2;end

if abs(out(2))<sqrt(eps) % mostlikely a line, p(1) is not reliable
    out(2) = 0; % to prevent scaling error in fc images
end
warning on

%%
% outlier rejection based on parametric model
yest = feval(model,out,x);
x_range = min(x):max(x);
y_range = feval(model,out,x_range);
outliers = abs(y-yest)>2;
%%
% if percentage of otliers is large, dont do correction!!
if sum(outliers)/length(outliers) < .25
    X_ = X_(~outliers,:);
    Y_ = Y_(~outliers,:);
    valid=1;
else
    out = out*0;
    valid = 0;
end

%%
if viz
    if iadj==1
        figure(304),
        subplot(2,1,iadj)
        cla
        plot(y,x,'+')
        hold on
        plot(y(outliers),x(outliers),'ro')
        plot(y_range,x_range,'g-')
        plot(feval(model,out,x),x,'go')
        daspect([1 100 1])
        
    elseif iadj==2
        figure(304),
        subplot(2,1,iadj)
        cla
        plot(x,y,'+')
        hold on
        plot(x(outliers),y(outliers),'ro')
        plot(x_range,y_range,'g-')
        plot(x,yest,'go')
        %     daspect([1 100 1])
    end
end
end

