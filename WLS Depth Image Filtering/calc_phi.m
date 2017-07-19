clc; clear; close all;
load('T0020_points.mat');

thetaH = pi/2;
thetaPrim = thetaH - theta ;
phiPrim = atan((z./y).*cos(thetaPrim));
Err = abs(phi - phiPrim);
imagesc(Err);figure(gcf);