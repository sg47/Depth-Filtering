function [inliers zPlane] = testFitRansac(X,Y,Z)
xx = X(:)';
yy = Y(:)';
zz = Z(:)';

[B P inliers] = ransacfitplane([xx; yy; zz],0.0004,0);
zPlane = - (B(4)+B(1)*xx(inliers) + B(2)*yy(inliers))/B(3);

xInliers = xx(inliers);
yInliers = yy(inliers);
zInliers = zz(inliers);
plot3(xInliers,yInliers,zPlane,'.b', 'markersize', 1);
hold on;
plot3(xInliers,yInliers,zInliers,'.r', 'markersize', 1);
return