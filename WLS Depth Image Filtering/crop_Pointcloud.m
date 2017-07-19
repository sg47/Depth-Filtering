function [ crop_x,crop_y,crop_z ] = crop_Pointcloud( x,y,z , x_min, x_max, y_min,y_max, z_min, z_max )
[m,n]=size(x);
crop_x=x; crop_y=y; crop_z=z;
for i=1:m
    for j=1:n
            if (~(x(i,j)>=x_min && x(i,j)<=x_max && y(i,j)>=y_min && y(i,j)<=y_max && z(i,j)>=z_min && z(i,j)<=z_max))
            %crop_x=[crop_x , x(i,j)];
            %crop_y=[crop_y , y(i,j)];
            %crop_z=[crop_z , z(i,j)];
            crop_x(i,j)=0; crop_y(i,j)=0; crop_z(i,j)=0;
            end
        
        end
end
end

