close all
clearvars
%% single frequency, power scan. plot current in couplers

%%config


freq = [8.5]*1e9;

X = [ 1.0709    1.0232    0.9961    1.0622    0.9548];

nockit_params = get_nockit6_params(X);
nockit_params.gnd_cond = 0;
input_idx= 4;
sig_pwr= -40; %dBm % for nockit6 params, critical power is ~-55 for couplers, and -32 for main lines
% iterations1 = 2*logspace(0,1,length(sig_pwr));
% iterations2 = 3*logspace(0,2,length(sig_pwr));

% iterations(1:(length(sig_pwr)-6)) = iterations1(1:(length(sig_pwr)-6));
% iterations((length(sig_pwr)-5):length(sig_pwr)) = iterations2((length(sig_pwr)-5):length(sig_pwr));
% iterations
iterations = [50]*ones(size(sig_pwr));
% critical_pwr = 10*log10((derived_params.Ic)^2/derived_params.Y0/1e-3);
% critical_pwr_c = 10*log10((derived_params.Icc)^2/derived_params.Y0/1e-3);


%%
%power loop:
% trans = zeros(length(sig_pwr),nockit_params.M, length(freq));
tic
for pwr_idx = 1:length(sig_pwr)
    [G, derived_params] = get_nockit_graph_NL(nockit_params, input_idx,sig_pwr(pwr_idx));   
    graph_data = process_graph_NL(G);
    
    
    txt_str = sprintf("calculating pwr #%d: %g dBm...", pwr_idx, sig_pwr(pwr_idx));
    disp(txt_str);

    [t_edges, r_edges] = solve_graph_NL_envelope(graph_data,freq, iterations(pwr_idx), true);
coordinates = get_nockit_coordinates(G,nockit_params);

% figure(501)
% % plt_network_power(G,coordinates, t_edges,r_edges);
% figure(502)
% plt_nockit_current(G,coordinates, t_edges,r_edges,freq, false,sig_pwr(pwr_idx));
% figure(503)
% plt_nockit_voltage(G,coordinates, t_edges,r_edges,freq,false,sig_pwr(pwr_idx));

figure(504)
% tiledlayout(1,1, 'padding', 'normal' )
plt_nockit_power(G,coordinates, t_edges,r_edges,freq,false, sig_pwr(pwr_idx));
end
disp("pwr scan finished")
toc


%%
