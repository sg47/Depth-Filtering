function ang = cosPatchNormal_4Side(nrm)
x = nrm(:,:,1);
y = nrm(:,:,2);
z = nrm(:,:,3);
[N_scan N_point] = size(x);
row_zero = zeros(1, N_point);
col_zero = zeros(N_scan, 1);
vertical(:,:,1) = [x(2:end,:); row_zero];
vertical(:,:,2) = [y(2:end,:); row_zero];
vertical(:,:,3) = [z(2:end,:); row_zero];
vertical_back(:,:,1) = [row_zero; x(1:end-1,:)];
vertical_back(:,:,2) = [row_zero; y(1:end-1,:)];
vertical_back(:,:,3) = [row_zero; z(1:end-1,:)];
horizontal(:,:,1) = [x(:, 2:end), col_zero];
horizontal(:,:,2) = [y(:, 2:end), col_zero];
horizontal(:,:,3) = [z(:, 2:end), col_zero];
horizontal_back(:,:,1) = [col_zero, x(:, 1:end-1)];
horizontal_back(:,:,2) = [col_zero, y(:, 1:end-1)];
horizontal_back(:,:,3) = [col_zero, z(:, 1:end-1)];

%% calc norm v and h of 4 side
out1 = dot(nrm,horizontal,3);
out1([1 N_scan], [1 N_point]) = NaN; 
out2 = dot(nrm,vertical,3);
out2([1 N_scan], [1 N_point]) = NaN; 
out3 = dot(nrm,horizontal_back,3);
out3([1 N_scan], [1 N_point]) = NaN;
out4 = dot(nrm,vertical_back,3);
out4([1 N_scan], [1 N_point]) = NaN;

ang=real(acos((out1+out2+out3+out4)/4));

end