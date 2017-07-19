function [cartScreen,sphScreen,angleStep] = sph(xyz, showflag)
% [m, n , ~] = size(xyz);
dim = 128;
i = 193:193+dim-1; j = 113:113+dim-1; m = dim; n = dim;

% X = xyz(:,:,1)'; % X = xyz(i,j,1)';%
% Y = xyz(:,:,2)'; % Y = xyz(i,j,2)';%
% Z = xyz(:,:,3)'; % Z = xyz(i,j,3)';%
 X = xyz(i,j,1)';%
 Y = xyz(i,j,2)';%
 Z = xyz(i,j,3)';%
x = X;
y = Z;
z = -Y;
[theta,phi,rho] = cart2sph(x,y,z);
dteta = abs(theta(n/2,m/2+1) - theta(n/2,m/2)); % horizontal space between angles
dphi = abs(phi(n/2+1,m/2) - phi(n/2,m/2));
angleStep = [dteta, dphi];% vertical space between angles
if showflag == 1
    plot3(x(:),z(:),y(:),'.b','markersize',1);
    xlabel('x');
    ylabel('y');
    zlabel('z');
end
cartScreen = cat(3,x,y,z);
sphScreen = cat(3,theta,phi,rho);