% Projects the given depth image to world coordinates. Note that this 3D
% coordinate space is defined by a horizontal plane made from the X and Z
% axes and the Y axis points up.
%
% Args:
%   imgDepthAbs - 480x640 depth image whose values indicate depth in
%                 meters.
%
% Returns:
%   X,Y,Z.
function [X Y Z] = depth_plane2depth_world(imgDepthAbs)
  [H, W] = size(imgDepthAbs);
 % assert(H == 480);
 % assert(W == 640);

  camera_params;

  [xx,yy] = meshgrid(1:W, 1:H);
  
  X = (xx ) .* imgDepthAbs / fx_d;
  Y = (yy) .* imgDepthAbs / fy_d;
  Z = imgDepthAbs;
   X = X/1000; Y = Y/1000; Z = Z/1000;
end
