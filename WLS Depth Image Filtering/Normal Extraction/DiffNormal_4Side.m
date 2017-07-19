function [ang] = DiffNormal_4Side(nrm)

[h, w, tmp] = size(nrm);
row_zero = zeros(1, w, 3);
col_zero = zeros(h, 1, 3);
vertical = [nrm(2:end, :, :); row_zero];
vertical_back = [row_zero; nrm(1:end-1, :, :)];
horizontal = [nrm(:, 2:end, :), col_zero];
horizontal_back = [col_zero, nrm(:, 1:end-1, :)];

%%  calc 4 normal vector
d_2x = dot(horizontal, horizontal_back, 3);
d_2y = dot(vertical, vertical_back, 3);

d_x = dot(horizontal, nrm, 3);
d_xb = dot(horizontal_back, nrm, 3);
d_y = dot(vertical, nrm, 3);
d_yb = dot(horizontal_back, nrm, 3);
ang = real(acos((2*d_2x+2*d_2y+d_x+d_xb+d_y+d_yb)/8));
tmp = max(max(ang(2:end-1,2:end-1)));
ang([1,h],:) = tmp; ang(:, [1,w]) = tmp;

end