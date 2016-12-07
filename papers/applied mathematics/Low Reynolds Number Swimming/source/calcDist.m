function dist = calcDist(matfile)
%function calcDist
eval(['load ' matfile]);

initialPt = zeros(3,41);

initialPt(1,:) = 0.4;
initialPt(2,:) = 0.2;
initialPt(3,:) = 0.4;

xVec(:,:) = xrods(1,:,:);
dif = xVec - initialPt;
size(dif)
dist(1:41) = sqrt(dif(1,:).^2 + dif(2,:).^2 + dif(3,:).^2);
