close all; clear all; clc

%1 flagellum
dist1 = calcDist1('newRK1');
save dist dist1;

close all; clear all; clc
%2 flagella
dist2 = calcDist1('newRK2');
load dist;
save dist dist1 dist2;

close all; clear all; clc
%4 flagella
dist4 = calcDist1('newRK4');
load dist;
save dist dist1 dist2 dist4;

close all; clear all; clc
%8 flagella
dist8 = calcDist1('newRK8');
load dist;
save dist dist1 dist2 dist4 dist8;

close all; clear all; clc

%now plot distances
time   = 0:.01:.4;
numPer = 1/(2*pi)*1/2*.1*time.^2;

load dist;
plot(numPer,dist1,'b',numPer,dist2,'r',numPer,dist4,'g',numPer,dist8,'m')
legend('1','2','4','8')
xlabel('Number of Revolutions')
ylabel('Distance Traveled (dim)');
title('Comparison Between Multiple Flagellar Configurations');