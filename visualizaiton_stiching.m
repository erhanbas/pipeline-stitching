matfolder = '/nrs/mouselight/cluster/classifierOutputs/2018-07-02/matfiles/'
matfolder = '/data3/visualization_stitching'
load(fullfile(matfolder,'regpts.mat'))
load(fullfile(matfolder,'scopeparams_pertile'),'paireddescriptor', ...
    'scopeparams', 'curvemodel','params')
load(fullfile(matfolder,'scopeloc'),'scopeloc','neighbors','experimentfolder','inputfolder')

addpath(genpath('./common'))
checkthese = [1 4 5 7]; % 0 - right - bottom - below
neigs = neighbors(:,checkthese);
% save(fullfile(matfolder,'vizdata.mat'))
targetfold_ = 'visualization_figures';

%%
clear ds
numtiles = length(paireddescriptor{1});
ds = zeros(numtiles,4);
for i=1:numtiles
    [ds(i,1:2)] = paireddescriptor{1}{i}.count;
    ds(i,3) = size(regpts{i}.X,1);
    if regpts{i}.matchrate>0
        ds(i,4) = regpts{i}.uni;
    end
end

%%
inds=(find(ds(:,1)>100&ds(:,2)>100&ds(:,3)>500))
[ds(inds,:) inds(:)]
val = ds(:,1).*ds(:,2);
[aa,bb] = find(val==max(val));
ineig = aa;
%%
model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model
dims = [1024 1536 251];
idxcent = ineig;
final_out = zeros(2,3);
xy=[];

for iadj = 1:2;
idxadj =  neigs(ineig,iadj+1);

if iadj==1
    descent = paireddescriptor{1}{ineig}.onx.X+1;
    descadj = paireddescriptor{1}{ineig}.onx.Y+1;
else
    descent = paireddescriptor{1}{ineig}.ony.X+1;
    descadj = paireddescriptor{1}{ineig}.ony.Y+1;
end


descent = util.correctTiles(descent,dims); % flip dimensions
descadj = util.correctTiles(descadj,dims); % flip dimensions
stgshift = 1000*(scopeloc.loc(idxadj,:)-scopeloc.loc(idxcent,:));
pixshift = round(stgshift.*(dims-1)./(params.imsize_um));
% descadj = descadj + ones(size(descadj,1),1)*pixshift; % shift with initial guess based on stage coordinate

outliers = zeros(size(descent,1),1);
iter = 1;
while true
    X_ = descent(~outliers,:);
    Y_ = descadj(~outliers,:);
    % figure,
    % myplot3(X_,'o')
    % hold on
    % myplot3(Y_,'o')
    
    dispvec = X_-Y_;
    y = dispvec(:,iadj);
    validinds = 1:length(y);
    X_ = X_(validinds,:);
    Y_ = Y_(validinds,:);
    y = y(validinds,:);
    
    x = X_(:,setdiff([1 2],iadj));
    dimcent = dims(setdiff([1:2],iadj))/2; % center of image along curvature axis
    
    out = curvemodel{1}(iadj,:,idxcent);
    options = optimoptions('fmincon','Display', 'off');
    options.ConstraintTolerance = 1e-9;
    options.OptimalityTolerance = 1e-9;
    fun = @(p) sum((y-feval(model,p,x)).^2);
    xyarr{iadj}.x = x;
    xyarr{iadj}.y = y;
    
    out2 = fmincon(fun,out,[],[],[],[],[],[],[],options);
    yest = feval(model,out2,sort(x));
    outliers = abs(y-yest)>=1;
    if mean(outliers)<.1
        break
    end
    iter =iter+1;
end
final_out(iadj,:) = out2;
end
%%
for iadj=1:2
    x = xyarr{iadj}.x;
    y = xyarr{iadj}.y;
if iadj == 1
    % close all
    figure(305),
    clf
    plot(y,x,'o','MarkerFaceColor','b')
    hold on
    % plot(x(outliers),y(outliers),'ro')
    % plot(x_range,y_range,'g-')
    plot(feval(model,final_out(iadj,:),-10:1500),-10:1500,'r-','LineWidth',4)
    axis equal
    daspect([1 10 1])
    xlim(final_out(iadj,3)+[-10 10])
    ylim([-10 1500])
    xlabel('\Delta X')
    ylabel('Y')
    export_fig(fullfile(targetfold_,'deltaX.png'),'-transparent')
    
else
    % close all
    figure(306),
    clf
    plot(x,y,'o','MarkerFaceColor','b')
    hold on
    % plot(x(outliers),y(outliers),'ro')
    % plot(x_range,y_range,'g-')
    plot(-10:1100,feval(model,final_out(iadj,:),-10:1100),'r-','LineWidth',4)
    axis equal
    daspect([10 1 1])

    ylim(final_out(iadj,3)+[-10 10])
    xlim([-10 1100])

    xlabel('\Delta Y')
    ylabel('Y')
    export_fig(fullfile(targetfold_,'deltaY.png'),'-transparent')
end
end
%%
order = params.order;
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

numTiles = length(scopeloc.filepath);
% xy crop for minimal overlap
[st,ed] = util.getcontolpixlocations(scopeloc,params,scopeparams);
[corrctrlpnttmp_,xlim_cntrl,ylim_cntrl] = util.getCtrlPts(params.imagesize(1),params.imagesize(2),params,st,ed);
% zlimdefaults = [5 25 dims(3)-26 dims(3)-6];
subcorrctrlpnttmp = corrctrlpnttmp_+1;
%%
[locs,xshift2D,yshift2D] = util.fcshift(final_out(:,:),order,xy,dims,[1 1]);
figure(34), 
imagesc(xshift2D)
caxis([-4 4])
axis equal tight
colorbar
    export_fig(fullfile(targetfold_,'xy_X.png'),'-transparent')


figure(43)
imagesc(yshift2D)
caxis([-2 2])
axis equal tight
colorbar
    export_fig(fullfile(targetfold_,'xy_Y.png'),'-transparent')

%%
% nbound = [0 0];
% nbound(1) = max(pixshift(iadj),min(descadj(:,iadj)));
% nbound(2) = min(dims(iadj),max(descent(:,iadj)))+0;
% [X_,Y_] = match.descriptorMatch(X,Y,params.matchparams);
paireddescriptor{1}{idxcent}.onx.X



%%
neigs_ = neighbors(ineig,:); %[2538 2539 2560 2797]
iter = 1;
clear Images
for iii  = neigs_([1 4 5 7])
    [pt,fi] = fileparts(scopeloc.filepath{iii});
    Images{iter} = mytiffread(fullfile(pt,sprintf('%s.0.tif',fi)));
    iter = iter + 1;
end

Ic = max(Images{1},[],3);
Ix = max(Images{2},[],3);
Iy = max(Images{3},[],3);
Iz = max(Images{4},[],3);

%%
% save stitching_vix Images ineig neigs_
%%
close all
figure(100),
imshow(Ic,[])
hold on
myplot3(paireddescriptor{2}{ineig}.onx.X+1,'bo'),
myplot3(paireddescriptor{2}{ineig}.ony.X+1,'bo'),
% myplot3(regpts{ineig}.X+1,'bo'),
% export_fig(fullfile(targetfold_,'center.png'))

%%
figure(101),
imshow(Ix,[])
hold on
myplot3(paireddescriptor{2}{ineig}.onx.Y+1,'ro'),
export_fig(fullfile(targetfold_,'Ix.png'))

%%
figure(102),
imshow(Iy,[])
hold on
myplot3(paireddescriptor{2}{ineig}.ony.Y,'ro'),
export_fig(fullfile(targetfold_,'Iy.png'))

%%
figure
imshow(Iz,[])
hold on,
myplot3(regpts{ineig}.Y,'r+')












%%
import // [73579.3, 19407.5, 31032.4]
%%
A = eye(4);
A(1:3,4)=[71067090,15072176,28216712]/1e3;
A(1:3,1:3) = diag([19384.43849206349,19238.116666666665,64652.963541666664])/2^6/1e3
%%
for i=1:2e4
    [ds(i,1:2)] = paireddescriptor{2}{i}.count;
end
%%
[bb]=find([ds(:,1).*ds(:,2)]==aa)
%%
ds(bb,:)

%%
ineig = 2538;
figure,
myplot3(paireddescriptor{2}{ineig}.onx.X,'bo'),
hold on,
myplot3(paireddescriptor{2}{ineig}.onx.Y,'r+')

% paireddescriptor{end}{tile_case}.onx
%%
checkthese = [1 4 5 7]; % 0 - right - bottom - below
imsize_um = params.imsize_um;
neigs =  neighbors(:,checkthese);%[id -x -y +x +y -z +z] format
%%

%%
x = paireddescriptor{2}{ineig}.onx.X;
y = paireddescriptor{2}{ineig}.onx.Y;
vec = x-y;

model = @(p,y) p(3) - p(2).*((y-p(1)).^2); % FC model

out = curvemodel{end}(1,:,ineig);

figure(304),
% subplot(2,1,iadj)
cla
plot(y,x,'+')
hold on
% plot(y(outliers),x(outliers),'ro')
% plot(y_range,x_range,'g-')
plot(feval(model,out,x),x,'go')
% daspect([1 100 1])

%%
for ii=1:1e5
    if scopeloc.filepath{ii} == '/groups/mousebrainmicro/mousebrainmicro/data/acquisition/2018-07-02/2018-07-06/01/01484/01484-ngc.acquisition'
        ii
        break
    end
end
