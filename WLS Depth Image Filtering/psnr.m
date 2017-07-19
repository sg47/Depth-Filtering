function psnr(rho, s, dh, dv, dd, sphScreen, angleStep,level, EPS,ext_mode)
%%
% hist(abs(dh{3}(:)),100);
% hist(abs(dv{3}(:)),100);
% hist(abs(dd{3}(:)),100);
%%
lvl = 10;
[m,n] = size(rho);
MSE = zeros(1,lvl);
PSNR = zeros(1,lvl);
MAX = max(rho(:));

% for lvl = level:-1:1
%     recExt_rho = iRDT2(s, dh, dv, dd, angleStep(1), angleStep(2), EPS, lvl, 4, 0,ext_mode);
%     rho_ = recExt_rho(1:m,1:n);
%     err = abs(rho-rho_);
%     err2 = err.^2;
%     MSE(lvl) = mean(err2(:));
%     PSNR(lvl) = 20*log10(MAX) - 10* log10(MSE(lvl));
%
%     [X,Y,Z] = sph2cart(sphScreen(:,:,1),sphScreen(:,:,2),rho_);
%     subplot(2,3,lvl);
%     plot3(X,Y,Z,'.b', 'markersize', 1);
%     title(int2str(lvl));
%
%     s = s(2:end);
%     dh = dh(2:end);
%     dv = dv(2:end);
%     dd = dd(2:end);
% end

for i= 1:lvl
    thr = i*0.0005;
    recExt_rho = iRDT2(s, dh, dv, dd, angleStep(1), angleStep(2), EPS, level, 4, thr,ext_mode);
    rho_ = recExt_rho(1:m,1:n);
    err = abs(rho-rho_);
    err2 = err.^2;
    MSE(i) = mean(err2(:));
    PSNR(i) = 20*log10(MAX) - 10* log10(MSE(i));

%     [X,Y,Z] = sph2cart(sphScreen(:,:,1),sphScreen(:,:,2),rho_);
%     subplot(2,3,i);
%     plot3(X,Y,Z,'.b', 'markersize', 1);
%     title(num2str(thr));
end
figure;
% subplot(211); plot(PSNR); title('PSNR');
% subplot(212); 
thr = 0.0005:0.0005:0.005;
plot(thr, MSE);title('MSE'); xlabel('Threshold value');

