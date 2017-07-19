clear;clc;
smooth_ = true;%false;%
close all;
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

Ray_XL = 555; Ray_XH = 570;
Ray_YL = 220; Ray_YH = 220;

Patch_XL = 550; Patch_XH = 580;
Patch_YL = 200; Patch_YH = 230;

%xx=xx(300:410,60:90); yy=yy(300:410,60:90); zz=zz(300:410,60:90);
%xx=xx(428:441,196); yy=yy(428:441,196); zz=zz(428:441,196); % plot index
% xx=xx(410:460,190:240); yy=yy(410:460,190:240); zz=zz(410:460,190:240);%2
% xx=xx(:); yy=yy(:); zz=zz(:);
% xx=xx(325:338); yy=yy(325:338); zz=zz(325:338); % plot index

%xx=xx(420:460,190:220); yy=yy(420:460,190:220); zz=zz(420:460,190:220);
%xx=xx(540:580,200:240); yy=yy(540:580,200:240); zz=zz(540:580,200:240); %3
%xx=xx(555:570,220); yy=yy(555:570,220); zz=zz(555:570,220); % ray patch3

s_nrm_4s = PatchNormal(xx,yy,zz);
s_nrm_4s=(1+s_nrm_4s)/2;
data(1,:)= xx(:); data(2,:)= yy(:); data(3,:)= zz(:);
color = ones(3,(size(xx,1)*size(xx,2))).*255;

% color_patch_x=ones(size(xx)); color_patch_y=ones(size(yy)); color_patch_z=ones(size(zz));
% color_patch_z(540:580,200:240)=0; color_patch_y(540:580,200:240)=0;
% color_patch(1,:,:)= color_patch_x(:); color_patch(2,:,:)= color_patch_y(:); color_patch(3,:,:)= color_patch_z(:);
% color_patch=round(mat2gray((color_patch)).*255);

%s_nrm_4s = reshape(s_nrm_4s, [3,size(s_nrm_4s,1)*size(s_nrm_4s,2)]);
nx=s_nrm_4s(:,:,1); ny=s_nrm_4s(:,:,2); nz=s_nrm_4s(:,:,3);
s_nrm(1,:,:)= nx(:); s_nrm(2,:,:)= ny(:); s_nrm(3,:,:)= nz(:);
s_nrm_4s = nonan(round(mat2gray((s_nrm)).*255));
%write_ply('data/sim_flat_xyz.ply',double(data),double(color));

x = xx; y = yy; z = zz;
[theta phi rho] = cart2sph(x,y,z);
img1=rho/max_range;
perfect_data = img1 * max_range;
[perfect_x perfect_y perfect_z] = sph2cart(theta,phi,img1* max_range);

RMS_out_array =[]; PSNR_out_array =[];
RMS_noise_array = []; PSNR_noise_array =[];
RMS_smooth_array = []; PSNR_smooth_array = [];
RMS_Bilateral_array = []; PSNR_Bilateral_array = [];
RMS_Combine_array = []; PSNR_Combine_array = [];

mean_array = [0.5];
variance_array = [0.00000002:0.000001:0.0001];
%variance_array = [0.00000002];
for i=1: 1 %size(variance_array,2)
    
        i
        distort = randn(size(x))*range_sigma;
        noise_img1 = img1 + distort;    % add random noise
        rho = rho + distort*max_range;
        
        v = variance_array(i)*var(img1(:));
        noise_img1 = imnoise(img1, 'gaussian', 0, v);
        %rho = rho *max_range;
        rho_max = max(rho(:));
        rho = rho / rho_max;
        rho = imnoise(rho, 'gaussian', 0, variance_array(i)).*rho_max;

        
        %save input_noise noise_img1 rho
        %load input_noise.mat;
        [x y z] = sph2cart(theta, phi, rho);
        %%


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
        %% normal edge regions detection 
        EdgeNormal=DiffNormal_4Side(s_nrm_4s);
        EdgeNormal(EdgeNormal>0.5) = 0.5;
        %%
        offsetVertical_Patch = (Ray_YL - Patch_YL )*(Patch_YH - Patch_YL + 1) + (Ray_XL  - Patch_XL +1);% 325;
        Ray_VerticalL=(Ray_YL-1)*size(x,1) + Ray_XL;
        Ray_VerticalH=(Ray_YH-1)*size(x,1) + Ray_XH;
        pointNum = Ray_VerticalH - Ray_VerticalL;
        
        
        Ray_VerticalL = Ray_VerticalL + 6;
        Ray_VerticalH = Ray_VerticalH - 7;
        offsetVertical_Patch =  offsetVertical_Patch + 7;
        pointNum = pointNum - 13;

        %% Trilateral filtering
        sigma = [space_sigma 10*range_sigma]; % bilateral filter standard deviations
        %tf_img1 = tfilter3(noise_img1,w,sigma,s_nrm_4s);

        bf_img1 = bfilter2(noise_img1(Patch_XL:Patch_XH,Patch_YL:Patch_YH),w,sigma);
        gf_img1 = gmdistFilter(noise_img1(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)),1);
        %gf_img1 = noise_img1(Patch_XL:Patch_XH,Patch_YL:Patch_YH);%gmdistFilter(noise_img1(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)),1);

        gfN_x = gmdistFilter(x(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)),0);
        gfN_y = gmdistFilter(y(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)),0);
        gfN_z = gmdistFilter(z(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)),0);
        %gL_img1 = gm_Laplacian(noise_img1(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)));

        gLN_x = gm_Laplacian_normal(x(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)));
        gLN_y = gm_Laplacian_normal(y(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)));
        gLN_z = gm_Laplacian_normal(z(Patch_XL:Patch_XH,Patch_YL:Patch_YH),nonan_Image(s_nrm_4s(Patch_XL:Patch_XH,Patch_YL:Patch_YH,:)));

        %gf_img1 = gL_img1;
        %% Mesh & Surf Plot

        %[tf_x, tf_y, tf_z] = sph2cart(theta, phi, tf_img1*max_range);
        [bf_x, bf_y, bf_z] = sph2cart(theta(Patch_XL:Patch_XH,Patch_YL:Patch_YH), phi(Patch_XL:Patch_XH,Patch_YL:Patch_YH), bf_img1*max_range);
        [gf_x, gf_y, gf_z] = sph2cart(theta(Patch_XL:Patch_XH,Patch_YL:Patch_YH), phi(Patch_XL:Patch_XH,Patch_YL:Patch_YH), gf_img1*max_range);
        %[theta_Out phi_Out rho_Out] = sph2cart(gLN_x, gLN_y, gLN_z);
        %[theta_Out phi_Out rho_Out] = sph2cart(gLN_x, gLN_y, gLN_z);
        
        EdgePatch = EdgeNormal(Patch_XL:Patch_XH,Patch_YL:Patch_YH);
        [x_Comb, y_Comb, z_Comb] = combineOutputs(gLN_x, gLN_y, gLN_z , gf_x, gf_y, gf_z , EdgePatch);
        
        
        
        
%         figure;
%         hold on;
% 
%         plot(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),'black','linewidth',2); 
%         %plot(x(Ray_VerticalL:Ray_VerticalH),y(Ray_VerticalL:Ray_VerticalH),'-black','linewidth',3,'markersize',10); 
%         plot(gfN_x(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),gfN_y(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),'-y','linewidth',1,'markersize',1); 
%         plot(gf_x(offsetVertical_Patch:offsetVertical_Patch+pointNum),gf_y(offsetVertical_Patch:offsetVertical_Patch+pointNum),'-r','linewidth',1,'markersize',5); 
%         plot(gLN_x(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),gLN_y(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),'-g','linewidth',1,'markersize',1); 
%         
%         %plot(tf_x(140828:140841),tf_y(140828:140841),'-y','linewidth',2);
%         %plot(bf_x(Ray_VerticalL:Ray_VerticalH),bf_y(Ray_VerticalL:Ray_VerticalH),'-g','linewidth',2);
%         legend('Noisy Image','gaussian normal','gaussian Filter - Laplacian(Proposed)' , 'Location','NorthWest');
%         xlabel('x(meter)');
%         ylabel('y(meter)');
%         print('-djpeg');
%         
        
        

        [perfect_theta perfect_phi perfect_rho] = sph2cart(perfect_x(Patch_XL:Patch_XH,Patch_YL:Patch_YH), perfect_y(Patch_XL:Patch_XH,Patch_YL:Patch_YH), perfect_z(Patch_XL:Patch_XH,Patch_YL:Patch_YH));
        [noise_theta noise_phi noise_rho] = sph2cart(x(Patch_XL:Patch_XH,Patch_YL:Patch_YH), y(Patch_XL:Patch_XH,Patch_YL:Patch_YH), z(Patch_XL:Patch_XH,Patch_YL:Patch_YH));
        [smooth_theta smooth_phi smooth_rho] = sph2cart(gf_x, gf_y, gf_z);
        [smoothNormal_theta smoothNormal_phi smoothNormal_rho] = sph2cart(gfN_x, gfN_y, gfN_z);
        [theta_Bilateral phi_Bilateral rho_Bilateral] = sph2cart(bf_x, bf_y, bf_z);

       
%         RMS_Out=calc_RMS(perfect_rho(offsetVertical_Patch:offsetVertical_Patch+pointNum),rho_Out(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1));
%         RMS_noise=calc_RMS(perfect_rho(offsetVertical_Patch:offsetVertical_Patch+pointNum),noise_rho(offsetVertical_Patch:offsetVertical_Patch+pointNum));
%         RMS_smooth=calc_RMS(perfect_rho(offsetVertical_Patch:offsetVertical_Patch+pointNum),smooth_rho(offsetVertical_Patch:offsetVertical_Patch+pointNum));
%       
        RMS_Out=calc_RMS_x_y(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),gLN_x(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),gLN_y(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1));
        RMS_noise=calc_RMS_x_y(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),x(Ray_VerticalL:Ray_VerticalH),y(Ray_VerticalL:Ray_VerticalH));
        RMS_smooth=calc_RMS_x_y(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),gf_x(offsetVertical_Patch:offsetVertical_Patch+pointNum),gf_y(offsetVertical_Patch:offsetVertical_Patch+pointNum));
        RMS_Combine=calc_RMS_x_y(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),x_Comb(offsetVertical_Patch:offsetVertical_Patch+pointNum),y_Comb(offsetVertical_Patch:offsetVertical_Patch+pointNum));

        %figure; plot (RMS_Out); hold on ; plot (RMS_smooth);

%         PSNR_Out=PSNR_Image(perfect_rho,rho_Out);
%         PSNR_noise=PSNR_Image(perfect_rho,noise_rho);
%         PSNR_smooth=PSNR_Image(perfect_rho,smooth_rho);


        PSNR_Out=PSNR_Image_xy(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),gLN_x(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1),gLN_y(offsetVertical_Patch-1:offsetVertical_Patch+pointNum-1));
        PSNR_noise=PSNR_Image_xy(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),x(Ray_VerticalL:Ray_VerticalH),y(Ray_VerticalL:Ray_VerticalH));
        PSNR_smooth=PSNR_Image_xy(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),bf_x(offsetVertical_Patch:offsetVertical_Patch+pointNum),bf_y(offsetVertical_Patch:offsetVertical_Patch+pointNum));
        PSNR_Combine=PSNR_Image_xy(perfect_x(Ray_VerticalL:Ray_VerticalH),perfect_y(Ray_VerticalL:Ray_VerticalH),x_Comb(offsetVertical_Patch:offsetVertical_Patch+pointNum),y_Comb(offsetVertical_Patch:offsetVertical_Patch+pointNum));

        RMS_out_array = [RMS_out_array RMS_Out];  RMS_noise_array = [RMS_noise_array RMS_noise];  RMS_smooth_array = [RMS_smooth_array RMS_smooth]; 
        PSNR_out_array = [PSNR_out_array PSNR_Out];  PSNR_noise_array = [PSNR_noise_array PSNR_noise];  PSNR_smooth_array = [PSNR_smooth_array PSNR_smooth]; 

end
%save arrays RMS_out_array RMS_noise_array RMS_smooth_array PSNR_out_array PSNR_noise_array PSNR_smooth_array variance_array;
plot_array(RMS_out_array,RMS_noise_array,RMS_smooth_array,PSNR_out_array,PSNR_noise_array,PSNR_smooth_array , variance_array);
% PSNR_noise=calc_RMS(perfect_z(Patch_XL:Patch_XH,Patch_YL:Patch_YH),z(Patch_XL:Patch_XH,Patch_YL:Patch_YH))
% PSNR_Out=calc_RMS(perfect_z(Patch_XL:Patch_XH,Patch_YL:Patch_YH),rho_Out)
%[theta phi rho] = cart2sph(x,z,-y);
%%



