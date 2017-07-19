function [normal cen] = PatchNormal_4Side(x, y, z)
[N_scan N_point] = size(x);
row_zero = zeros(1, N_point);
col_zero = zeros(N_scan, 1);
vertical(:,:,1) = x - [x(2:end,:); row_zero];
vertical(:,:,2) = y - [y(2:end,:); row_zero];
vertical(:,:,3) = z - [z(2:end,:); row_zero];
vertical_back(:,:,1) = x - [row_zero; x(1:end-1,:)];
vertical_back(:,:,2) = y - [row_zero; y(1:end-1,:)];
vertical_back(:,:,3) = z - [row_zero; z(1:end-1,:)];
horizontal(:,:,1) = x - [x(:, 2:end), col_zero];
horizontal(:,:,2) = y - [y(:, 2:end), col_zero];
horizontal(:,:,3) = z - [z(:, 2:end), col_zero];
horizontal_back(:,:,1) = x - [col_zero, x(:, 1:end-1)];
horizontal_back(:,:,2) = y - [col_zero, y(:, 1:end-1)];
horizontal_back(:,:,3) = z - [col_zero, z(:, 1:end-1)];

%%  calc 4 normal vector
normal1 = cross(horizontal, vertical, 3);
normal2 = cross(vertical_back, horizontal, 3);
normal3 = cross(horizontal_back, vertical_back, 3);
normal4 = cross(vertical, horizontal_back, 3);

% normal1(N_scan, N_point, :) = 1; 
% normal2(1, N_point, :) = 1; 
% normal3(1, 1, :) = 1; 
% normal4(N_scan, 1, :) = 1;
normal1([1 N_scan], [1 N_point], :) = 1;
normal2([1 N_scan], [1 N_point], :) = 1;
normal3([1 N_scan], [1 N_point], :) = 1;
normal4([1 N_scan], [1 N_point], :) = 1;
normal1 = normal1./repmat(sqrt(sum(normal1.^2, 3)), [1 1 3]);                 % Length normalization
normal2 = normal2./repmat(sqrt(sum(normal2.^2, 3)), [1 1 3]);                 % Length normalization
normal3 = normal3./repmat(sqrt(sum(normal3.^2, 3)), [1 1 3]);                 % Length normalization
normal4 = normal4./repmat(sqrt(sum(normal4.^2, 3)), [1 1 3]);                 % Length normalization

normal = (normal1+normal2+normal3+normal4)/4;
normal = normal./repmat(sqrt(sum(normal.^2, 3)), [1 1 3]);                 % Length normalization

cen = zeros(size(normal));
cen(:,:,1) = (5*x - horizontal(:,:,1) - vertical(:,:,1) - vertical_back(:,:,1) - horizontal_back(:,:,1))/5;
cen(:,:,2) = (5*y - horizontal(:,:,2) - vertical(:,:,2) - vertical_back(:,:,2) - horizontal_back(:,:,2))/5;
cen(:,:,3) = (5*z - horizontal(:,:,3) - vertical(:,:,3) - vertical_back(:,:,3) - horizontal_back(:,:,3))/5;

normal = repmat(-sign(dot(normal,cen,3)),[1 1 3]).*normal;
end