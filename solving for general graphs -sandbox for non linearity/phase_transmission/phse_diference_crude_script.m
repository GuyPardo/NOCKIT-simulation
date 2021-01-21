
clearvars phase

for line = 1:4
freq_scan_nockit6_phase_transmission11
phase(line,:) = trans_phase(line,:);




end

phase_delta = phase(1,:) - phase(4,:);
%%

figure(706)
clf
for line=1:4
subplot(2,2,line)
plot(freq, phase(line,:))
P = polyfit(1e-9*freq, phase(line,:),1 );
yfit = 1e-9*P(1)*freq +P(2);
slope = 1e-9*P(1);
total_L = L0*(N-1) + 2*lengths(line);
v_ph_est = 2*pi*total_L/slope;
 
hold on
plot(freq, yfit, '--')
title(sprintf("line %g , v_p_h_,_e_s_t = %.4g",line, v_ph_est))
grid on
end

figure(705)
clf
    plot(freq, phase_delta, 'linewidth', 1.5)
    grid on
    leg = legend(num2str((line)'),"location", "best", "fontsize", 16);
    
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("phase", "fontsize", 16)
    title_str = sprintf('line 1 - line 4  \n NOCKIT5 fit = %d \n coplanar couplers = %d' ,  nockit5_fit,coplanar_couplers);
    title(title_str, "fontsize", 16)

    
    P = polyfit(1e-9*freq, phase_delta,1 );
    yfit = 1e-9*P(1)*freq +P(2);
    slope = 1e-9*P(1);
    
    total_L = 2*(lengths(1) - lengths(4));
    v_ph_est = 2*pi*total_L/slope;
    
    hold on
    plot(freq, yfit, '--')
    leg = legend([sprintf('1-->1 - 4-->4' ), sprintf("linear fit, slope = %.3g\n v_p_h_,_e_s_t = %.4g", slope, v_ph_est)],"location", "best", "fontsize", 16);
    
