function [ PSNR_Value ] = PSNR_Image_xy( Actual_x , Actual_y , x , y )

[rows columns ~] = size(Actual_x);

% Calculate mean square error of R, G, B.   
err = sqrt((Actual_x-x).^2 +(Actual_y-y).^2) ;

% Then "square" the "error".
squareError = err.^2;

mseRImage = squareError;
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

