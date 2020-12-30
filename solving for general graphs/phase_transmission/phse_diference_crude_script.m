line = 1;
freq_scan_nockit6_phase_transmission11
phase1 = trans_phase(line,:);
line = 4;
freq_scan_nockit6_phase_transmission11
phase4 = trans_phase(line,:);

phase_delta = phase1-phase4;


figure(705)
clf
    plot(freq, phase_delta, 'linewidth', 1.5)
    grid on
    leg = legend(num2str((line)'),"location", "best", "fontsize", 16);
    
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("phase", "fontsize", 16)
    title_str = sprintf('phase at output w. input from line %d \n NOCKIT5 fit = %d \n coplanar couplers = %d' , input_idx, nockit5_fit,coplanar_couplers);
    title(title_str, "fontsize", 16)

    
    P = polyfit(1e-9*freq, phase_delta,1 );
    yfit = 1e-9*P(1)*freq +P(2);
    slope = 1e-9*P(1);
    
    total_L = 2*(lengths(1) - lengths(4));
    v_ph_est = 2*pi*total_L/slope;
    
    hold on
    plot(freq, yfit, '--')
    leg = legend([sprintf('1-->1 - 4-->4' ), sprintf("linear fit, slope = %.3g\n v_p_h_,_e_s_t = %.4g", slope, v_ph_est)],"location", "best", "fontsize", 16);
    
