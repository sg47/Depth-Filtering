clc;
close all;
img1 =double(rho)/max(rho(:));
figure
imagesc(img1);
title('input image');
%%
 ret = bilateral_me (img1, 5, 1, 1,C);
 figure
%imagesc(ret-img1);
imagesc(ret);
title('result of bilateral filter')
 %%
