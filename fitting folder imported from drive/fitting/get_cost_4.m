function [cost] = get_cost_4(freq, data_dB, x)

% measurement data:
tic

    

     
display(x)
    
     %default  config
%     N = 31;
%     L = 100e-6;
%     d = 27e-6;
%     M = 7;
    

    % default config
%     v_ph  = 1.6850e+06;
%     v_ph_c = 1.7626e+06;
%     Z_0 = 49.8391;
%     Z_c = 547.9262;
    
    
    %updatd default config
    
%     v_ph = 1.3611e+06;
%     v_ph_c = 1.4088e+06;
%     Z_0 = 49.3596;
%     Z_c = 365.6212;
    
    
    v_ph = 1.8628e+06;
    v_ph_c = 5.4916e+06;
    Z0 = 50.1153;
    Zc = 1.9916e+03;
    
    
    v_ph_new = v_ph*x(1);
    v_ph_c_new = v_ph_c*x(2);
    Zc_new = Zc*x(3);
 
    
   trans = freq_scan_fun_4(freq,v_ph_new,v_ph_c_new,1/Z0,1/Zc_new );
   trans_dB = 20*log10(abs(trans));
   %size(trans)
   %size(data_dB)
    cost = sum(sum(abs(trans_dB(4,:) + x(4) - data_dB(4,:)).^2));

toc
end
