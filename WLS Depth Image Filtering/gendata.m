
clc
clear;

width = 128;
height = 128;

[x z] = meshgrid(1:width, 1:height);
x = x - width/2;
z = z - height/2;
z = -z;

x = x/10;
z = z/10;

y = 0*x + 0*z + 4;

rgb = ones(height,width,3)*128;

[theta phi rho] = cart2sph(x,y,z);

save('plane.mat', 'x', 'y', 'z', 'theta', 'rho', 'phi');