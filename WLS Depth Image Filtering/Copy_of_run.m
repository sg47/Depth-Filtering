clear;clc;
smooth_ = true;%false;%

addpath 'data';
addpath 'depth2point';
addpath 'Normal Extraction';
addpath 'normalSimilarity';
addpath 'ply convertor';

%%
max_range = 25;
w = 10;
space_sigma = 5;
range_sigma = 0.02/max_range; % 0.02/max_range;

%% Loading data

load sim_flat_xyz.mat;
xx = x-4; zz = y-1.7; yy = -z-5;            % to center the sensor

%xx=xx(300:410,60:90); yy=yy(300:410,60:90); zz=zz(300:410,60:90);
%xx=xx(428:441,196); yy=yy(428:441,196); zz=zz(428:441,196); % plot index
% xx=xx(410:460,190:240); yy=yy(410:460,190:240); zz=zz(410:460,190:240);
% xx=xx(:); yy=yy(:); zz=zz(:);
% xx=xx(325:338); yy=yy(325:338); zz=zz(325:338); % plot index

%xx=xx(420:460,190:220); yy=yy(420:460,190:220); zz=zz(420:460,190:220);

s_nrm_4s = PatchNormal(xx,yy,zz);
s_nrm_4s=(1+s_nrm_4s)/2;
data(1,:)= xx(:); data(2,:)= yy(:); data(3,:)= zz(:);
color = ones(3,(size(xx,1)*size(xx,2))).*255;

%s_nrm_4s = reshape(s_nrm_4s, [3,size(s_nrm_4s,1)*size(s_nrm_4s,2)]);
nx=s_nrm_4s(:,:,1); ny=s_nrm_4s(:,:,2); nz=s_nrm_4s(:,:,3);
s_nrm(1,:,:)= nx(:); s_nrm(2,:,:)= ny(:); s_nrm(3,:,:)= nz(:);
s_nrm_4s = nonan(round(mat2gray((s_nrm)).*255));
%write_ply('data/sim_flat_xyz.ply',double(data),double(color));

x = xx; y = yy; z = zz;
[theta phi rho] = cart2sph(x,y,z);
img1=rho/max_range;

distort = randn(size(x))*range_sigma;

noise_img1 = img1 + distort;    % add random noise
rho = rho + distort*max_range;
[x y z] = sph2cart(theta, phi, rho);
%%

if smooth_ == true
    %h = fspecial('gaussian',10,15);
    %h = fspecial('gaussian',10,4);
    h = fspecial('gaussian',10,25);
    s_rho = imfilter(rho,h);
    [s_x, s_y, s_z] = sph2cart(theta, phi, s_rho);
    s_nrm = PatchNormal(s_x,s_y,s_z);
    s_nrm_4s = PatchNormal_4Side(s_x,s_y,s_z);
%    s_angle_normal=cosPatchNormal(s_nrm);
%    s_angle_normal_4Side=cosPatchNormal_4Side(s_nrm_4s);
    
    
 %   figure(1); subplot(1,2,1); imshow((1+s_nrm)/2);title('normal')
 %   subplot(1,2,2); imshow((1+s_nrm_4s)/2);title('4side-normal')
%    subplot(1,4,3); imshow(mat2gray(s_angle_normal)); title('Cos angle')
%    subplot(1,4,4); imshow(mat2gray(s_angle_normal_4Side)); title('4side-Cos angle')
else
    nrm = PatchNormal(x,y,z);
    nrm_4s = PatchNormal_4Side(x,y,z);
%    angle_normal=cosPatchNormal(nrm);
%    angle_normal_4Side=cosPatchNormal_4Side(nrm_4s);
    figure(2); subplot(1,2,1); imshow((1+nrm)/2);title('normal')
    subplot(1,2,2); imshow((1+nrm_4s)/2);title('4side-normal')
%    subplot(1,4,3); imshow(mat2gray(angle_normal)); title('Cos angle')
%    subplot(1,4,4); imshow(mat2gray(angle_normal_4Side)); title('4side-Cos angle')
end


%% Trilateral filtering
sigma = [space_sigma 10*range_sigma]; % bilateral filter standard deviations
%tf_img1 = tfilter3(noise_img1,w,sigma,s_nrm_4s);
bf_img1 = bfilter2(noise_img1,w,sigma);

gm_img1 = gmdistFilter(noise_img1(410:460,190:240),nonan_Image(s_nrm_4s(410:460,190:240,:)));


%% Mesh & Surf Plot

%[tf_x, tf_y, tf_z] = sph2cart(theta, phi, tf_img1*max_range);
[bf_x, bf_y, bf_z] = sph2cart(theta, phi, bf_img1*max_range);
[gf_x, gf_y, gf_z] = sph2cart(theta(410:460,190:240), phi(410:460,190:240), gm_img1*max_range);

xxx=x; yyy=y; zzz=z;
bf_xxx=bf_x; bf_yyy=bf_y; bf_zzz=bf_z;
%tf_xxx=tf_x; tf_yyy=tf_y; tf_zzz=tf_z;
[xx yy zz]= crop_Pointcloud(xx,yy,zz,-8.8,-7,-8,-5.4,-1,2);
[gf_x gf_y gf_z]= crop_Pointcloud(gf_x,gf_y,gf_z,-8.8,-7,-8,-5.4,-1,2);

[x y z]= crop_Pointcloud(xxx,yyy,zzz,-8.8,-7,-8,-5.4,-1,2);
[bf_x bf_y bf_z]= crop_Pointcloud(bf_xxx,bf_yyy,bf_zzz,-8.8,-7,-8,-5.4,-1,2);
%[tf_x tf_y tf_z]= crop_Pointcloud(tf_xxx,tf_yyy,tf_zzz,-8.8,-7,-8,-5.4,-1,2);


figure;
hold on;
offset = 325;
pointNum = 13;
plot(xx(140828:140841),yy(140828:140841),'black','linewidth',2); 
plot(x(140828:140841),y(140828:140841),'-*b','linewidth',2,'markersize',10); 
plot(gf_x(offset:offset+pointNum),gf_y(offset:offset+pointNum),'-*r','linewidth',2,'markersize',10); 

%plot(tf_x(140828:140841),tf_y(140828:140841),'--y','linewidth',2);
plot(bf_x(140828:140841),bf_y(140828:140841),'-go','linewidth',2);
legend('Orginal Image','Noisy Image','guassian normal(proposed)','Bilateral Filter');
xlabel('x(meter)');
ylabel('y(meter)');
print('-djpeg');
%[theta phi rho] = cart2sph(x,z,-y);
%%



