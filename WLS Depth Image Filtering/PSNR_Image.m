function [ PSNR_Value ] = PSNR_Image( I,Ihat )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
% Read the dimensions of the image.
[rows columns ~] = size(I);

% Calculate mean square error of R, G, B.   
mseRImage = (double(I(:,:,1)) - double(Ihat(:,:,1))) .^ 2;
%mseGImage = (double(I(:,:,2)) - double(Ihat(:,:,2))) .^ 2;
%mseBImage = (double(I(:,:,3)) - double(Ihat(:,:,3))) .^ 2;
for i=1:rows
    for j=1:columns
        if(isnan(mseRImage(i,j))) mseRImage(i,j)=0;
        end
    end
end
mseR = sum(sum(mseRImage)) / (rows * columns);
%mseG = sum(sum(mseGImage)) / (rows * columns);
%mseB = sum(sum(mseBImage)) / (rows * columns);

% Average mean square error of R, G, B.
%mse = (mseR + mseG + mseB)/3;

% Calculate PSNR (Peak Signal to noise ratio).
PSNR_Value = 10 * log10( 255^2 / mseR);

end

