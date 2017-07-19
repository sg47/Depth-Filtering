function ang = cosPatchNormal(nrm)
x = nrm(:,:,1);
y = nrm(:,:,2);
z = nrm(:,:,3);
[N_scan N_point] = size(x);
row_zero = zeros(1, N_point);
col_zero = zeros(N_scan, 1);
vertical(:,:,1) = [x(2:end,:); row_zero];
vertical(:,:,2) = [y(2:end,:); row_zero];
vertical(:,:,3) = [z(2:end,:); row_zero];
horizontal(:,:,1) = [x(:, 2:end), col_zero];
horizontal(:,:,2) = [y(:, 2:end), col_zero];
horizontal(:,:,3) = [z(:, 2:end), col_zero];

out1 = dot(nrm,horizontal,3);
out1([1 N_scan], [1 N_point]) = NaN; 
out2 = dot(nrm,vertical,3);
out2([1 N_scan], [1 N_point]) = NaN; 

ang=real(acos((out1+out2)/2));

end