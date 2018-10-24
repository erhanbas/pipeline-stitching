function edgelist = drawPlates(xyz,edgelist,input_col,grid_format)
% creates a plate from given control points
if nargin<4
    grid_format = [5 5 4]; % [numof control points in X and Y and Z]
end
if size(xyz,1) ~=prod(grid_format)
    error('check grid size')
end

numOfControlPointsinLayer = (grid_format(1)*grid_format(2));
if nargin<2 | isempty(edgelist)
    for ix =1:grid_format(3)
        inds = [1:numOfControlPointsinLayer] + (ix-1)*numOfControlPointsinLayer;
        xyz_layer = xyz(inds,:);
        for ixx = 1:2:grid_format(1)-1
            for idx_i = 1:floor(grid_format(1)/2)
                idx_control = grid_format(1)*ixx+2*idx_i;
                xyz_layer(idx_control,1) = xyz_layer(idx_control,1) + .1;
                xyz_layer(idx_control+1,1) = xyz_layer(idx_control+1,1) - .1;
            end
        end
        DT = delaunayTriangulation(xyz_layer(:,1:2));
        T = triangulation(DT.ConnectivityList,xyz_layer);
        edgelist{ix} = T.ConnectivityList + grid_format(1)*grid_format(2)*(ix-1);
    end
    edgelist = cat(1,edgelist{:});
end
%%
% add vertical edges
[xyz_layer1,xyz_layer2]=deal([]);
for ix =1:grid_format(3)-1
    inds1 = [1:numOfControlPointsinLayer] + (ix-1)*numOfControlPointsinLayer;
    inds2 = [1:numOfControlPointsinLayer] + (ix)*numOfControlPointsinLayer;
    xyz_layer1{ix} = xyz(inds1,:);
    xyz_layer2{ix} = xyz(inds2,:);
end
xyz_layer1 = cat(1,xyz_layer1{:});
xyz_layer2 = cat(1,xyz_layer2{:});
X = [xyz_layer1(:,1) xyz_layer2(:,1) NaN(size(xyz_layer1(:,1)))];
Y = [xyz_layer1(:,2) xyz_layer2(:,2) NaN(size(xyz_layer1(:,2)))];
Z = [xyz_layer1(:,3) xyz_layer2(:,3) NaN(size(xyz_layer1(:,3)))];

if nargout>0
    return
end
myfig = figure(gcf); 
clc
mygca = gca;
set(myfig,'Color',[1 1 1]*.1)
set(mygca,'Color',[1 1 1]*.5)

% create gradient on given color 
c = mean(reshape(xyz(edgelist,3),[],3),2);
c = (c-min(c))/(max(c)-min(c));
[grad,im]=colorGradient(input_col,input_col*.5,256);
% interpolate grad using c
col = interp1(linspace(1,0,256),grad,c);
h = patch('faces',edgelist,'vertices',xyz,'FaceVertexCData',col,...
      'facecolor','flat', ...
      'edgecolor',[0 0 0],'parent',mygca,...
      'FaceAlpha',.8);

% trisurf(edgelist,xyz(:,1),xyz(:,2),xyz(:,3),'edgecolor',[0 0 0],'FaceAlpha',0.8)
hold on
plot3(X',Y',Z','color','k')
set(gca,'Zdir','reverse')
% axis equal 
% axis tight


