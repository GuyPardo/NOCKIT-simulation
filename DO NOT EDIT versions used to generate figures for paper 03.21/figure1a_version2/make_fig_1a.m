
clearvars
load NOCKIT5_2traces_data.mat;
freq_red = downsample(freq,401);
% X = [t,W,Wc,H,lam]

X = [0.8341    0.9311    1.2711     1.0139    0.6633]; % latest fit result as of 30.5.21

get_nockit2_params;

G = change_params(nockit_params,X);

offset = -21; %dB


graph_data = process_graph(G);

trans = nan(2,length(freq_red));
% loop on frequencies
for i=1:length(freq_red)
    [t_edges, r_edges] = solve_graph(graph_data,freq_red(i)); % solve
    
  
    
    [trans(:,i), ~] = read_nockit_solution(nockit_params, G, t_edges, r_edges);
    
end
%%fff
   trans_dB = 20*log10(abs(trans));
   trans_mag2 = abs(trans);
   data_mag2 = abs(data).*(10^(+2.1/2));
%%
figure(607); clf;
  xlabel('Frequency (MHz)','fontsize', 15); ylabel(sprintf('Normalized\nTransmission [a.u.]'), 'fontsize', 15); 
 hold on
 cc = colororder();

  for i=1:2
   
    plot(1e-6*freq, data_mag2(:,i), 'color',cc(i,:), 'linewidth', 2)
  end

for i=1:2
 plot(1e-6*freq_red, (trans_mag2(i,:)), 'linewidth',2 ,'linestyle', '--','color',cc(i,:));
end
 
 legend(["Measured 1 \rightarrow 2", "Measured 1 \rightarrow 4", "Simulated 1 \rightarrow 2","Simulated 1 \rightarrow 4" ], 'location', 'southwest', 'fontsize', 15)
set(gca, 'YScale', 'log')
ylim([1e-4, 2])
% grid on
