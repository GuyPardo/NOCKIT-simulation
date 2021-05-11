clearvars
folder = 'Z:\Measurements\2020.dir\01_12_20 cooldown\Nockit6\Transmission@-40\';

for i=1:4
    ss = sprintf("4to%d",i); 
    path = strcat(folder, ss);
    [freq(i,:), data(i,:)] = import_data_4_fun(path);
    
   
end

freq = freq(1,:);

data_dB = 20*log10(abs(data));
%%
 figure(808)
plot(freq,20*log10(abs(data)), 'linewidth',2); grid on; title ("measured transmission w. input from line 4"); ylabel("dB", 'fontsize',16);xlabel("frequency (Hz)", 'fontsize',16);
leg = legend(num2str((1:4)'),"location", "best", "fontsize", 13);
title(leg, "line")

%%
smooth_num=30;
for i=1:4
     data_smooth_mag(i,:) = smooth((abs(data(i,:))),smooth_num );
end

data_smooth_dB = 20*log10(data_smooth_mag);

figure(809)
plot(freq,20*log10(data_smooth_mag), 'linewidth',2); grid on; title ("smoothed transmission w. input from line 4"); ylabel("dB", 'fontsize',16);xlabel("frequency (Hz)", 'fontsize',16);
leg = legend(num2str((1:4)'),"location", "best", "fontsize", 13);
title(leg, "line")

