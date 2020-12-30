
% clearvars


coplanar_couplers = true;
nockit5_fit  = false;
%% construct graph
% this part constructs the graph representing nockit network: a square
% lattice made out of  M lines with N couplers between each adjacent pair. 



tic
% geometry: and network structure
N=31; % number of couplers. (= number of unit cells minus 1) 
M = 7; % number of lines
L0 = 100e-6; % length of each line segment
d = 27e-6; % length of each coupler segment
t = 10e-9;
W = 3e-6;
W_c = 200e-9;
H = 16e-9;
gap_c = 1.4e-6;
addpath(genpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation'))

[Y0, v_ph]  = get_microstrip_properties(W,t,H);
if coplanar_couplers
    [Yc, v_ph_c ] = get_CPW_properties(t,W_c,gap_c);
else
    [Yc, v_ph_c ] = get_microstrip_properties(W_c, t,H);
end

% % % physical parameters: (see the NOCKIT simulation for the way these values were calculated  )
% % v_ph = 1.361104539023962e+06; % phase velocity for lines
% % v_ph_c =  1.408763793738406e+06; % phase velocity for couplers
% % Y0 = 0.020259488114653; % admittance for lines
% % Yc = 0.002735070881675; % admittance for couplers

if nockit5_fit
% parameters correction from fit. use these to get somthing close to the
% measurement for 2 traces NOCKIT5, but note that we still have to explain the factor of 2 in
% the phase velocity. the other two factors are close to 1, so they are OK.
x = [1.9935    0.9193    0.8418];     
%x = [2,1,1]
v_ph = v_ph*x(1);
    v_ph_c = v_ph_c*x(2);
    Yc =  Yc/x(3);
    
end


% frequency etc.
freq = 1e9*linspace(3,9,201); 
omega= 2*pi*freq;



input_idx = [line];   % can be more than one.
%%
% define graph: define an array of nodes with M rows and N+2 columns. the
% nodes are numbered such that nodes 1:M are the first column, M+1:2*M are
% the socond column  etc.
nodes = reshape(1:M*(N+2),M,N+2 );  
G = digraph();
% define main lines:
% connect each node (m,n) to (m, n+1), with weight  2.
% by convention weight 2 means a "main" edge and not a coupler.
s = nodes(:,1:N+1);
t = nodes(:,2:N+2);
w = 2*ones(size(s));
G = G.addedge(s,t,w);

% define coupling lines:
% connect each node (m,n) to (n+1,n) with weight 1.
% by convention weight 1 means a coupling edge.
s = reshape(nodes(1:M-1, 2:N+1), 1,(N)*(M-1));
t = reshape(nodes(2:M, 2:N+1), 1,(N)*(M-1));
w = ones(size(s));
G = G.addedge(s,t,w);


edge_num = G.numedges;
nodes_num = G.numnodes;

% define coordinates for plotting the graph: (this has no effect on the solution)
x = repmat(L0*(0:N+1), M,1);
y = repmat(d*fliplr((0:M-1)), 1,N+2); % the y coordinates are in the flipped to plot from top to bottom

x = reshape(x, 1,nodes_num);
y = reshape(y, 1,nodes_num);
LWidths = 6*G.Edges.Weight/max(G.Edges.Weight); 
figure(504)
clf
G.plot('xdata', x, 'ydata',y, 'linewidth', LWidths);

%


% define edges attributres; phase velocity, length and characteristic
% admittance:
% rememeber weight = 1 means coupler edge, weight = 2 means regular edge:
clearvars G.Edges.v_ph G.Edges.L G.Edges.Y
G.Edges.v_ph(G.Edges.Weight==2) = v_ph*ones(sum(G.Edges.Weight==2),1);
G.Edges.v_ph(G.Edges.Weight==1) = v_ph_c*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.Y(G.Edges.Weight==2) = Y0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Y(G.Edges.Weight==1) = Yc*ones(sum(G.Edges.Weight==1),1);

% define boundary conditions attribute for each edge according to the
% following convention:
% 0 - do nothing
% 1 - set t to 0
% 2 - set r to 0
% 3 - set t to 1
% 4 - set r t0 1
G.Edges.BC = zeros(edge_num,1);
G.Edges.BC(G.findedge(nodes(:,1),nodes(:,2))) = 1*ones(M,1);
G.Edges.BC(G.findedge(nodes(input_idx,1),nodes(input_idx,2))) = 3;
G.Edges.BC(G.findedge(nodes(:,N+1),nodes(:,N+2))) = 2*ones(M,1);


lengths = [2600,2000, 1700, 500, 1700,2000,2600]*1e-6; 




G.Edges.L(G.findedge(nodes(:,1),nodes(:,2))) =lengths'; 
G.Edges.L(G.findedge(nodes(:,end-1),nodes(:,end))) =lengths';






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
    ref(:,i) = r_edges(G.findedge(nodes(:,1),nodes(:,2)));
    trans(:,i) =t_edges(G.findedge(nodes(:,end-1),nodes(:,end))); 
    
    
    
end
toc
%%
ref_mag2 = abs(ref);
trans_mag2 = abs(trans);

 trans_phase = unwrap(angle(trans), [],2);

%% plot
ref_dB = 10*log10(ref_mag2);
trans_dB = 10*log10(trans_mag2);

figure(701) % transmittance graphs
clf;
    plot(freq, trans_dB(1:4,:), 'linewidth', 1.5)
    grid on
    leg = legend(num2str((1:4)'),"location", "best", "fontsize", 13);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("dB", "fontsize", 16)
    title("transmittace", "fontsize", 16)
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
 %%      
figure(704)
clf
    plot(freq, trans_phase(line,:), 'linewidth', 1.5)
    grid on
    leg = legend(num2str((line)'),"location", "best", "fontsize", 16);
    
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("phase", "fontsize", 16)
    title_str = sprintf('phase at output w. input from line %d \n NOCKIT5 fit = %d \n coplanar couplers = %d' , input_idx, nockit5_fit,coplanar_couplers);
    title(title_str, "fontsize", 16)

    
    P = polyfit(1e-9*freq, trans_phase(line,:),1 );
    yfit = 1e-9*P(1)*freq +P(2);
    slope = 1e-9*P(1);
    
    total_L = L0*(N-1) + 2*lengths(line);
    v_ph_est = 2*pi*total_L/slope;
    
    hold on
    plot(freq, yfit, '--')
    leg = legend([sprintf('%d-->%d', line,line), sprintf("linear fit, slope = %.3g\n v_p_h_,_e_s_t = %.4g", slope, v_ph_est)],"location", "best", "fontsize", 16);
    
