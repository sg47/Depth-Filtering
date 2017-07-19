function write_ply(fname, Depth, Color)
color = Color;
range = Depth;
%[x y z] = depth_plane2depth_world(double(range));
cameraMatrix = getNYU2Camera();
[x y z] = getPointCloudFromZ(double(range), cameraMatrix, 1);

[normal cen] = PatchNormal_4Side(x, y, z);
color = uint8(reshape(mat2gray(color),size(color,1)*size(color,2),3).* 255);
color = color';
data = zeros(3,(size(x,1)*size(x,2)));
%color = ones(3,(size(x,1)*size(x,2))).*255;
for i=1:size(x,1)
    for j=1:size(x,2)
        index = i + (j-1)*size(x,1);
        data(1,index) = x(i,j);
        data(2,index) = y(i,j);
        data(3,index) = z(i,j);
    end
end
write_ply(fname,data,color);
