function newMap = newColMap(centerPoint,scalingIntensity,dataMin,dataMax);
if nargin<2
    centerPoint = 25;
    scalingIntensity = 5;
    dataMin = 0;
    dataMax = 255;
end
cMap = redbluecmap(256);
x = 1:length(cMap);
x = x - (centerPoint-dataMin)*length(x)/(dataMax-dataMin);
x = scalingIntensity * x/max(abs(x));
x = sign(x).* exp(abs(x));
x = x - min(x); x = x*511/max(x)+1;
newMap = interp1(x, cMap, 1:512);