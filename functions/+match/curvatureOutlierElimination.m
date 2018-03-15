function [paireddescriptor,curvemodel,unreliable] = curvatureOutlierElimination(paireddescriptor,curvemodel,scopeloc)
%% outlier rejection based on median model
Nneig = size(curvemodel,3);
validtiles=squeeze(all(curvemodel(1:2,1,:)|curvemodel(1:2,3,:)));
[thrs,medcurvemodel] = estimateFCthreshols(curvemodel(:,:,validtiles),.99);
medcurvemodelcents = medcurvemodel(1:2,[1 3]);
thrcurvemodelcents = thrs(1:2,[1 3]);
% count = zeros(Nneig,2);
% for ineig = 1:Nneig
%     [count(ineig,:)] = paireddescriptor{ineig}.count;
% end
unreliable = zeros(Nneig,1);
for ineig = 1:Nneig
    fccent = curvemodel(1:2,[1 3],ineig);
    if any(abs(fccent(:)-medcurvemodelcents(:))>thrcurvemodelcents(:))
        unreliable(ineig) = 1;
    end
end

inliers = find(~unreliable);
% for every tiles estimate an affine
anchors = scopeloc.gridix(inliers,1:3);
queries = scopeloc.gridix(:,1:3);
IDX = knnsearch(anchors,queries,'K',1,'distance',@distfun);%W=[1 1 100000]

% fill missing 
for ineig = 1:Nneig
    ianch = inliers(IDX(ineig));
    paireddescriptor{ineig}.onx = paireddescriptor{ianch}.onx;
    paireddescriptor{ineig}.ony = paireddescriptor{ianch}.ony;
    paireddescriptor{ineig}.count = [size(paireddescriptor{ineig}.onx.X,1) size(paireddescriptor{ineig}.ony.X,1)];
    curvemodel(:,:,ineig) = curvemodel(:,:,ianch);
end
