function [ RMS ] = calc_RMS( Actual , Predicted )
err = Actual - Predicted;
 
% Then "square" the "error".
squareError = err.^2;
 
% Then take the "mean" of the "square-error".
meanSquareError = mean(mean(squareError));
 
% Then take the "root" of the "mean-square-error" to get 
% the root-mean-square-error!
 
RMS = sqrt(meanSquareError);

end

