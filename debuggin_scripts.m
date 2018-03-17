load(fullfile(matfolder,'scopeparams_pertile_singletile'),'paireddescriptor', ...
    'scopeparams', 'R', 'curvemodel','scopeparams_', 'paireddescriptor_', ...
    'curvemodel_','params')
%%
[X_,Y_] = match.descriptorMatch(X,Y,matchparams);
Y_(:,iadj) = Y_(:,iadj)- pixshift(iadj);% move it back to original location after CDP
% get field curvature model
[X_,Y_,out] = match.fcestimate(X_,Y_,iadj,matchparams);



%%
figure
myplot3(X_,'.')
hold on
myplot3(Y_,'.')

%%
icent = 5136
ineig = 5163

% icent = 5461
% ineig = 5488

Acent = -scopeparams_(icent).affineglFC/1e3;
Aneig = -scopeparams_(ineig).affineglFC/1e3;
Xcent = paireddescriptor_{icent}.ony.X;
Ycent = paireddescriptor_{icent}.ony.Y;

order = params.order;
dims = params.imagesize;
xlocs = 1:dims(1);
ylocs = 1:dims(2);
[xy2,xy1] = ndgrid(ylocs(:),xlocs(:));
xy = [xy1(:),xy2(:)];

[Xcent_fc] = util.fcshift(curvemodel(:,:,icent),order,xy,dims,Xcent+1);
Xcent_fc = Xcent_fc-1;

[Ycent_fc] = util.fcshift(curvemodel(:,:,icent),order,xy,dims,Ycent+1);
Ycent_fc = Ycent_fc-1;

shift = median(Xcent-Ycent);
shift_fc = median(Xcent_fc-Ycent_fc);
locX = [[Acent 1e3*scopeloc.loc(icent,:)']*[Xcent ones(size(Xcent,1),1)]']';
locY = [[Aneig 1e3*scopeloc.loc(ineig,:)']*[Ycent ones(size(Xcent,1),1)]']';
locX_fc = [[Acent 1e3*scopeloc.loc(icent,:)']*[Xcent_fc ones(size(Xcent,1),1)]']';
locY_fc = [[Aneig 1e3*scopeloc.loc(ineig,:)']*[Ycent_fc ones(size(Xcent,1),1)]']';

figure(icent)
dcm = datacursormode( gcf );
dcm.UpdateFcn = @NewCallback;
clf
subplot(311)
hold on
myplot3(locX,'b.')
myplot3(locY,'ro')
axis equal
title(sprintf('%d/%d/%d',shift))
subplot(312)
hold on
myplot3(locX_fc,'b.')
myplot3(locY_fc,'ro')
axis equal
title(sprintf('%d/%d/%d',round(shift_fc)))

figure(icent*10)
dcm = datacursormode( gcf );
dcm.UpdateFcn = @NewCallback;
clf
subplot(311)
hold on
myplot3(Xcent,'b.')
myplot3(Ycent+shift,'ro')
axis equal
title(sprintf('%d/%d/%d',shift))
subplot(312)
hold on
myplot3(Xcent_fc,'b.')
myplot3(Ycent_fc+shift_fc,'ro')
axis equal
title(sprintf('%d/%d/%d',round(shift_fc)))

%%
figure(icent)
dcm = datacursormode( gcf );
dcm.UpdateFcn = @NewCallback;
clf
subplot(311)
hold on
myplot3(Xcent,'b.')
myplot3(Ycent+shift,'ro')
axis equal
title(sprintf('%d/%d/%d',shift))
subplot(312)
hold on
myplot3(Xcent_fc,'b.')
myplot3(Ycent_fc+shift_fc,'ro')
axis equal
title(sprintf('%d/%d/%d',round(shift_fc)))
datacursormode(gcf)

%%



%%
% for debugging scripts
ineig = 5163
tifname = scopeloc.filepath{ineig}
filefold = fileparts(tifname);
[~,tifname] = fileparts(filefold);
tifpath = fullfile(filefold,sprintf('%s-ngc.0.tif',tifname));
Ic = mytiffread(tifpath);
%
try
    ineigx = neighbors(ineig,4)
    ineigy = neighbors(ineig,5)
catch
    ineigx = neigs(ineig,2)
    ineigy = neigs(ineig,3)
end

tifname = scopeloc.filepath{ineigy}
filefold = fileparts(tifname);
[~,tifname] = fileparts(filefold);
tifpath = fullfile(filefold,sprintf('%s-ngc.0.tif',tifname));
In = mytiffread(tifpath);

%%
ImC = max(Ic,[],3);
ImN = max(In,[],3);

ImX = rot90(ImC,2);
ImY = rot90(ImN,2);
%%
curvemodel(:,3,[5162 5163 5164])
%%
figure(7)
cla
imagesc(density,'Xdata',[gridx],'Ydata',[gridy])
axis tight
hold on,
plot(x,y,'m.')

%%
figure(11), cla
imshow(ImX,[11 15]*1e3)
hold on
myplot3(descent,'d')
myplot3(X,'o')
myplot3(X_,'go')
% myplot3(X_(prob_inliers,:),'go')

%%
figure(21), cla
imshow(ImY,[11 15]*1e3)
hold on
myplot3(descadj-pixshift,'d')
myplot3(Y-pixshift,'o')
myplot3(Y_,'+')

%%
dispvec = X_-Y_;
figure,
plot(X_(:,1),dispvec(:,2),'.')
% xlim([min(dispvec(:,1)) max(dispvec(:,1))])
axis tight

%%
figure(1), 
cla
imshow(ImC,[11 15]*1e3)
hold on
myplot3(paireddescriptor{ineig}.onx.X,'o')
% myplot3(descriptors{ineig}(:,1:3),'d')
% camroll(90)


%%
figure(2), 
imshow(ImN,[11 15]*1e3)
hold on
myplot3(paireddescriptor{ineig}.ony.Y,'o')

%%
figure(2), 
imshow(ImXp1,[])
hold on
myplot3(Y_,'o')
% myplot3(descent,'o')
%%
old=[-376.7212665974943206, -1.4150539397464144, 44.3289718325040667, 0.0, 71108753.59580696,
    2.3374766378038374, -297.0775556009028833, -13.8336760551395486, 0.0, 14667134.00371697,
    0.2542090768871866, 0.7845732492211573, 888.2112152288162861, 0.0, 35412268.12609629,
    0.0, 0.0, 0.0, 1.0, 0.0, 
    -0.0, 0.0, 0.0, 0.0, 1.0]
new=[-376.2621493364757725, -1.4195567316836457, 44.7492455332472261, 0.0, 71108071.94444494,
    2.2485697966564264, -296.7169950905342830, -14.2602659659113709, 0.0, 14667817.50149242,
    0.0439235442969926, 0.5983385231528975, 887.8041013352451500, 0.0, 35413245.98378776,
    0.0, 0.0, 0.0, 1.0, 0.0,
    -0.0, 0.0, 0.0, 0.0, 1.0];

round(old*[1;1;1;0;1]/1000)'
round(new*[1;1;1;0;1]/1000)'
round(old*[1024;1536;1;0;1]/1000)'
round(new*[1024;1536;1;0;1]/1000)'
%%
err=[]
for itt = 33:35
    er1=scopeloc.loc(find(scopeloc.gridix(:,1)==220&scopeloc.gridix(:,2)==itt&scopeloc.gridix(:,3)==779),:)
    err(end+1,:) = er1(1,:);
end
%%
err3=[]
for ix = 218:220
for iy = 33:35
    er1=scopeloc.loc(find(scopeloc.gridix(:,1)==ix&scopeloc.gridix(:,2)==iy&scopeloc.gridix(:,3)==781),:)
    err3(end+1,:) = er1(1,:);
end
end
diff(err3)
%%
er2=scopeloc.loc(find(scopeloc.gridix(:,2)==33&scopeloc.gridix(:,3)==780),:)
er3=scopeloc.loc(find(scopeloc.gridix(:,2)==34&scopeloc.gridix(:,3)==780),:)

