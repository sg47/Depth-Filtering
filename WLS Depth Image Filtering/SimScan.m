function [cartScreen,sphScreen,angleStep] = SimScan(anglesInfo,FOV, showflag)
% [cartScreen,sphScreen,angleStep] = SimScan([-60,60,-45,45],[256,256]);
tic;
leftAngle = anglesInfo(1);
rightAngle = anglesInfo(2);
downAngle = anglesInfo(3);
topAngle = anglesInfo(4);
row = FOV(1);
col = FOV(2);
hAngles = linspace(leftAngle,rightAngle,col); % horizontal angles of screen
vAngles = linspace(downAngle,topAngle,row); % vertical angles of screen
hAngles = (pi/180)* hAngles;
vAngles = (pi/180)* vAngles;
dteta = hAngles(2) - hAngles(1); % horizontal space between angles
dphi = vAngles(2) - vAngles(1);
angleStep = [dteta, dphi];% vertical space between angles
%%
% hAngles = pi/2 - hAngles;
minDistance = 10;
x = zeros(length(vAngles),length(hAngles));
z = x;
y = x;
for i = 1:length(vAngles)
        x(i,:) = minDistance*tan(hAngles)./cos(vAngles(i));%%%%% edit
        z(i,:) = -ones(1,length(hAngles)).*(minDistance*tan(vAngles(i)));
end
y(:) = 10;
%%
[theta,phi,rho]=cart2sph(x,y,z);
% noise = (rand(size(rho)) - 0.5)*0.05;
noise = 0;
% noise=zeros(size(rho));
% noise(row/2-2:row/2+2,col/2-2:col/2+2)=(randn(5)-.5)*.05;
% % rho = rho + noise;
% % [x, y, z] = sph2cart(theta,phi,rho);
%%
if showflag == 1
    plot3(x(:),z(:),y(:),'.b');
    xlabel('x');
    ylabel('y');
    zlabel('z');
end

cartScreen = cat(3,x,y,z);
sphScreen = cat(3,theta,phi,rho);
fprintf('generation time = %f\n',toc);

