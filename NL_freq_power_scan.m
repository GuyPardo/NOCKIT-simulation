%% config

%freq scan
freq = linspace(3,9,201)*1e9;
% get parameters

nockit_params = get_nockit6_params();
input_idx = 4;
sig_pwr = -80; % in dBm

G = get_nockit_graph_NL(nockit_params, input_idx, sig_pwr);

