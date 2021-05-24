%%config

freq = linspace(3,9,101)*1e9;

X = [ 1.0709    1.0232    0.9961    1.0622    0.9548];

nockit_params = get_nockit6_params(X);
input_idx= 1;
sig_pwr= linspace(-80,-35,20); %dBm % for nockit6 params, critical power is ~-55 for couplers, and -32 for main lines
iterations1 = 2*logspace(0,1,length(sig_pwr));
iterations2 = 3*logspace(0,2,length(sig_pwr));

iterations(1:(length(sig_pwr)-6)) = iterations1(1:(length(sig_pwr)-6));
iterations((length(sig_pwr)-5):length(sig_pwr)) = iterations2((length(sig_pwr)-5):length(sig_pwr));
% iterations

% critical_pwr = 10*log10((derived_params.Ic)^2/derived_params.Y0/1e-3);
% critical_pwr_c = 10*log10((derived_params.Icc)^2/derived_params.Y0/1e-3);


%%
%power loop:
trans = zeros(length(sig_pwr),nockit_params.M, length(freq));
tic
for pwr_idx = 1:length(sig_pwr)
    [G, derived_params] = get_nockit_graph_NL(nockit_params, input_idx,sig_pwr(pwr_idx));   
    graph_data = process_graph_NL(G);
    
    
    txt_str = sprintf("calculating pwr #%d: %g dBm...", pwr_idx, sig_pwr(pwr_idx));
    disp(txt_str);
    for i=1:length(freq)
       [t_edges, r_edges] = solve_graph_NL_envelope(graph_data,freq(i), iterations(pwr_idx));

       [t,r]   = read_nockit_solution(nockit_params, G,t_edges,r_edges);

       trans(pwr_idx,:,i) = t(:,end)./derived_params.sig_amp; % normalized transmission

    end

end
disp("pwr scan finished")
toc
%%
trans_dB = 20*log10(abs(trans));
trans_dB_norm = trans_dB - repmat(trans_dB(1,:,:), length(sig_pwr),1,1);
%% Plot


output_idx = 1; 
clear plt_mat
plt_mat(:,:) = trans_dB(:,output_idx,:);


figure(301)
clf
plot(freq,plt_mat, 'linewidth', 2); grid on

 cmap = colormap(jet(length(sig_pwr))) ; %Create Colormap
 colororder(cmap);
 cbh = colorbar ; %Create Colorbar
 N_ticks = 10;
 cbh.Ticks = linspace(0, 1,N_ticks) ; 
 cbh.TickLabels = num2cell(linspace(min(sig_pwr), max(sig_pwr),N_ticks)) ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
cbh.Label.String = "signal power (dBm)";
cbh.Label.FontSize = 16;
xlabel("frequency (Hz)", 'fontsize', 15)
ylabel("dB", 'fontsize', 15)
title(sprintf("transmission %d-->%d", input_idx, output_idx));


figure(302); clf; surf(plt_mat, 'Xdata', freq, "Ydata",sig_pwr )
view(2), shading flat
xlabel("frequency (Hz)", 'fontsize', 15)
ylabel("signal power (dBm)", 'fontsize', 15)
title(sprintf("transmission %d-->%d", input_idx, output_idx));
colormap default
cbh = colorbar;
cbh.Label.String = "dB";
cbh.Label.FontSize = 16;
axis('xy')

%
clear plt_mat
plt_mat(:,:) = trans_dB_norm(:,output_idx,:);


figure(303)
clf
plot(freq,plt_mat, 'linewidth', 2); grid on

 cmap = colormap(jet(length(sig_pwr))) ; %Create Colormap
 colororder(cmap);
 cbh = colorbar ; %Create Colorbar
 N_ticks = 10;
 cbh.Ticks = linspace(0, 1,N_ticks) ; 
 cbh.TickLabels = num2cell(linspace(min(sig_pwr), max(sig_pwr),N_ticks)) ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
cbh.Label.String = "signal power (dBm)";
cbh.Label.FontSize = 16;
xlabel("frequency (Hz)", 'fontsize', 15)
ylabel(sprintf("s_%d_%d - s_%d_%d(P_m_i_n)",output_idx,input_idx,output_idx,input_idx), 'fontsize', 15)
title(sprintf("normalized transmission %d-->%d", input_idx, output_idx));


figure(304); clf; surf(plt_mat, 'Xdata', freq, "Ydata",sig_pwr )
view(2), shading flat
xlabel("frequency (Hz)", 'fontsize', 15)
ylabel("signal power (dBm)", 'fontsize', 15)
title(sprintf("normalized transmission %d-->%d", input_idx, output_idx));
colormap default
cbh = colorbar;
cbh.Label.String = sprintf("s_%d_%d - s_%d_%d@P_m_i_n",output_idx,input_idx,output_idx,input_idx);
cbh.Label.FontSize = 16;
% axis('xy')


%%

saveplots = false;

if saveplots
    
    folder  = "C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\repos2\NOCKIT-simulation\NL solver\figs\";
    filename = sprintf("sim_%dto%d_1D_%d_pwrs_max_iter_%d", input_idx, output_idx, length(sig_pwr), iterations(end));
    figure(301); savefig(strcat(folder, filename))
    
    filename = sprintf("sim_%dto%d_2D_%d_pwrs_max_iter_%d", input_idx, output_idx, length(sig_pwr), iterations(end));
    figure(302); savefig(strcat(folder, filename))
    
    filename = sprintf("sim_%dto%d_norm_1D_%d_pwrs_max_iter_%d", input_idx, output_idx, length(sig_pwr), iterations(end));
    figure(303); savefig(strcat(folder, filename))
    
    filename = sprintf("sim_%dto%d_norm_2D_%d_pwrs_max_iter_%d", input_idx, output_idx, length(sig_pwr), iterations(end));
    figure(304); savefig(strcat(folder, filename))
end