function trans = freq_scan_non_lin_fun(freq, v_ph, v_ph_c,Y0, Yc, Ic, Icc, sig_pwr,d)
% here we create a nockit2 grpah and solving with the non linear solver
% 
%% construct graph
% this part constructs the graph representing the 2 traces ladder network
% (NOCKIT) 


gnd_cond = 0; % loss: conductance to ground;
% intereference experiment setup
%  N_pwrs = 20;
% sig_pwr =  -75; % dbm

iterations = 0;


% geometry: and network structure
N=30; % number of couplers. (= number of unit cells minus 1) 
M = 2; % number of lines
L0 = 100e-6; % length of each line segment
% d = 27e-6; % length of each coupler segment
thickness = 9e-9;%10e-9;
W = 2.3e-6;
W_c = 300e-9;
H = 16e-9;
gap_c = 1.4e-6;


input_idx = [1]; 


gnd_cond_c = gnd_cond*W_c/W;
addpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation')


C = Y0/v_ph;
Cc = Yc/v_ph_c;
L = 1/v_ph/Y0;
Lc = 1/v_ph_c/Yc;
% %Ic = 0.0002*t_new/0.00000001*W/0.000001; % samuel's formula
% Ic = 60e-6*(t_new/6e-9)*(W/2.3e-6); % % from mikita measurement
% % Ic= 1000000000; % if we only want to consider the coupler's non linearity
% Icc = 0.0002*t_new/0.00000001*W_c_new/0.000001; % samuel's formula
% Icc = 60e-6*(t_new/6e-9)*(W_c_new/2.3e-6); % from mikita measurement


sig_amp = sqrt(1/Y0*10^((sig_pwr/10) - 3)); % in Volts

% frequency etc.
%freq = linspace(freq_min,freq_max,freq_res); % in Hz 
omega= 2*pi*freq;
k0 = 2*pi*freq/v_ph; % wavenumber for lines
kc = 2*pi*freq/v_ph_c; % wavenumber for couplers

% define grpah: define an array of nodes with M rows and N+2 columns. the
% nodes are numbered such that nodes 1:M are the first column, M+1:2*M are
% the socind column  etc.
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


%

% define edges attributres; pahe velocity, length and characteristic
% admittance:
% rememeber weight = 1 means coupler edge, weight = 2 means regular edge:
clearvars G.Edges.v_ph G.Edges.L G.Edges.Y
G.Edges.C(G.Edges.Weight==2) = C*ones(sum(G.Edges.Weight==2),1);
G.Edges.C(G.Edges.Weight==1) = Cc*ones(sum(G.Edges.Weight==1),1);
G.Edges.Len(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Len(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = Lc*ones(sum(G.Edges.Weight==1),1);
G.Edges.Ic(G.Edges.Weight==2) = Ic*ones(sum(G.Edges.Weight==2),1);
G.Edges.Ic(G.Edges.Weight==1) = Icc*ones(sum(G.Edges.Weight==1),1);
G.Edges.gnd_cond(G.Edges.Weight==2) = gnd_cond*ones(sum(G.Edges.Weight==2),1);
G.Edges.gnd_cond(G.Edges.Weight==1) = gnd_cond_c*ones(sum(G.Edges.Weight==1),1);

% define boundary conditions attribute for each edge according to the
% following convention:
% 0 - do nothing
% 1 - set t to 0
% 2 - set r to 0
% 3 - set t to G.Edges.BCval
% 4 - set r t0 G.Edges.BCval
G.Edges.BC = zeros(edge_num,1);
G.Edges.BC(G.findedge(nodes(:,1),nodes(:,2))) = 1*ones(M,1);
G.Edges.BC(G.findedge(nodes(input_idx,1),nodes(input_idx,2))) = 3;
G.Edges.BC(G.findedge(nodes(:,N+1),nodes(:,N+2))) = 2*ones(M,1);

G.Edges.BCval = zeros(edge_num,1);
G.Edges.BCval(G.findedge(nodes(input_idx(1),1),nodes(input_idx(1),2))) = sig_amp;



% L_arr = G.Edges.L;
% C_arr = G.Edges.C;
% Len_arr = G.Edges.Len;
% BC_arr = G.Edges.BC;
% BCval_arr = G.Edges.BCval;
% %

graph_data = process_graph(G);

%% loop on freqs
trans = zeros(M,length(freq));
for i = 1:length(freq)
    [t_edges, r_edges] = solve_graph_non_lin_2(graph_data,freq(i), iterations);
    
    % read solution: this part is specific to the NOCKIT geometry 
 t = reshape(t_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
 r = reshape(r_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
    
    
    trans(:,i) = t(:,end)/sig_amp;
    
end
% trans_dB = 20*log10(abs(trans));

end
%%
