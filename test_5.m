
clearvars


% geometry: and network structure


nockit_params.N=31; % number of couplers. (= number of unit cells minus 1) 
nockit_params.M = 7; % number of lines
nockit_params.L0 = 100e-6; % length of each line segment
nockit_params.d = 27e-6; % length of each coupler segment
nockit_params.t = 8.5e-9; % thickness of metal
nockit_params.W = 3e-6; % width of main traces
nockit_params.W_c = 200e-9; % width of couplers
nockit_params.H = 35e-9; % thickness of dielectric
nockit_params.gap_c = 8e-6; % for coplanar couplers, width of the gap between trace and ground.
nockit_params.input_idx = [4];

G = get_nockit_graph(nockit_params);

M = nockit_params.M;
N = nockit_params.N;

freq = linspace(2,9,801)*1e9;

%%  solve
% pre-process grpah:
graph_data = process_graph(G);

% pre-allocate:
ref = nan(M,length(freq));
trans = nan(M,length(freq));
% loop on frequencies
for i=1:length(freq)
    [t_edges, r_edges] = solve_graph(graph_data,freq(i)); % solve
    
  % read solution: (this part is specific to the NOCKIT geometry)   
[ref(:,i),trans(:,i)] = read_nockit_solution(nockit_params, G, t_edges, r_edges);
    
end
toc

ref_mag2 = abs(ref);
trans_mag2 = abs(trans);

trans_phase = unwrap(angle(trans));

%% plot
ref_dB = 10*log10(ref_mag2);
trans_dB = 10*log10(trans_mag2);

figure(701) % transmittance graphs
clf;
    plot(freq, trans_dB(1:4,:), 'linewidth', 2)
    grid on
    leg = legend(num2str((1:4)'),"location", "best", "fontsize", 13);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("dB", "fontsize", 16)
    title(sprintf("transmittace w. input from line %d", nockit_params.input_idx), "fontsize", 16)
figure(702)  % reflectance graphs  
clf;
    plot(freq, ref_dB, 'linewidth', 1.5)
    grid on
    leg = legend(num2str((1:M)'),"location", "best", "fontsize", 13);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("dB", "fontsize", 16)
    title("reflectance", "fontsize", 16)

figure(703) % colormaps
clf;
   subplot(2,1,1) % transmittance
       imagesc( trans_dB, "XData", freq, "YData", (1:M))
       view(2)
       shading flat
       colorbar
       xlabel("frequency (Hz)", "fontsize", 13)
       ylabel("line", "fontsize", 16)
       title("transmittace (dB)", "fontsize", 16)  
    subplot(2,1,2) % reflectance
       imagesc( ref_dB, "XData", freq, "YData", (1:M))
       view(2)
       shading flat
       colorbar
       xlabel("frequency (Hz)", "fontsize", 13)
       ylabel("line", "fontsize", 16)
       title("reflectance (dB)", "fontsize", 16)
    
       colormap jet

       