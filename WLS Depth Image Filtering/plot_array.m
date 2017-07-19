function [ output_args ] = plot_array (RMS_out_array,RMS_noise_array,RMS_smooth_array,PSNR_out_array,PSNR_noise_array,PSNR_smooth_array , variances);

figure;
hold on;
plot(variances, RMS_out_array , 'r' );
plot(variances, RMS_noise_array , 'g' );
plot(variances,RMS_smooth_array , 'b' );
legend('RMS Normal filter','RMS Noise','RMS gaussian Filter' , 'Location','NorthWest');
xlabel('Variance value of noise');
ylabel('RMS');
hold off;
%% 
figure;
hold on;
plot(variances,PSNR_out_array , 'r' );
plot(variances, PSNR_noise_array , 'g' );
plot(variances,PSNR_smooth_array , 'b' );
legend('PSNR Normal filter','PSNR Noise','PSNR gaussian Filter' , 'Location','NorthWest');
xlabel('Variance value of noise');
ylabel('PSNR');
hold off;

end

