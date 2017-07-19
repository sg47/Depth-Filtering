function [x_Comb, y_Comb, z_Comb] = combineOutputs(gLN_x, gLN_y, gLN_z , bf_x, bf_y, bf_z , EdgePatch);

EdgePatch = mat2gray(EdgePatch);
[theta_N phi_N rho_N] = cart2sph(gLN_x, gLN_y, gLN_z);
[theta_B phi_B rho_B] = cart2sph(bf_x, bf_y, bf_z);

landa = 0.5;
%alpha = EdgePatch;
alpha = exp(-1 .* landa .* EdgePatch);

rho_out = alpha.* rho_N + (1-alpha).* rho_B;

[x_Comb, y_Comb, z_Comb] = sph2cart(theta_N, phi_N, rho_out);

end

