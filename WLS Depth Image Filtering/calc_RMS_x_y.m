function [ RMS ] = calc_RMS_x_y( Actual_x , Actual_y , x , y )
err = sqrt((Actual_x-x).^2 +(Actual_y-y).^2) ;

% Then "square" the "error".
squareError = err.^2;
 
% Then take the "mean" of the "square-error".
% meanSquareError = mean(mean(squareError));
meanSquareError = sum(sum(squareError));
 
% Then take the "root" of the "mean-square-error" to get 
% the root-mean-square-error!
 
RMS = sqrt(meanSquareError);
%RMS = err;
end
