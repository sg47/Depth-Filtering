clc
%clear all
%close all
%I=imread ('untitled.png');
%I=rgb2gray(I);
%I=I>130;
I=B;
regmax = imregionalmax(I,8);
imshow(regmax);
figure ; imagesc(I);
%% example 1 
% A = 10*ones(10,10);
% A(2:4,2:4) = 22; 
% A(6:8,6:8) = 33; 
% A(2,7) = 44;
% A(3,8) = 45;
% A(4,9) = 44;
% A
% regmax = imregionalmax(A)
