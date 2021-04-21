% sandbox
% here we create a nockit2 grpah and solving with the non linear solver
% 
clearvars
close all
 nockit5_fit=true;
 coplanar_couplers=false;
%% construct graph
% this part constructs the graph representing the 2 traces ladder network
% (NOCKIT) 


gnd_cond = 0; % loss: conductance to ground;
% intereference experiment setup
N_pwrs = 1;
sig_pwr =  -55; % dbm
pump_min = -32;
pump_max = -21;
pump_pwr = linspace(pump_min,pump_max,N_pwrs); % dbm

iterations = 1*logspace(1,2,N_pwrs);
% iterations = linspace(10,100,N_pwrs);
 iterations = 67*ones(size(pump_pwr));

phase = linspace(0,2*pi, 20);
figure(204)
clf
colororder(jet(N_pwrs));


% geometry: and network structure
N=31; % number of couplers. (= number of unit cells minus 1) 
M = 2; % number of lines
L0 = 100e-6; % length of each line segment
d = 27e-6; % length of each coupler segment
thickness = 8e-9;%10e-9;
W = 2.3e-6;
W_c = 300e-9;
H = 16e-9;
gap_c = 1.4e-6;


trans = zeros(N_pwrs, length(phase), M);
trans_norm = zeros(N_pwrs, length(phase), M);
for pwr_idx = 1:N_pwrs
    for ii = 1:length(phase)

sig_amp = sqrt(50*10^((sig_pwr/10) - 3))*exp(-1i*phase(ii)); % in Volts, assuming here it's about 50 ohm.\
pump_amp = sqrt(50*10.^((pump_pwr(pwr_idx)/10) - 3));


input_idx = [2,1]; % should be a vector of length two. the first component is the signal idx, the second is the pump

tic


gnd_cond_c = gnd_cond*W_c/W;
addpath(genpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation'))

[Y0, v_ph]  = get_microstrip_properties(W,thickness,H);
if coplanar_couplers
    [Yc, v_ph_c ] = get_CPW_properties(thickness,W_c,gap_c);
else
    [Yc, v_ph_c ] = get_microstrip_properties(W_c, thickness,H);
end

C = Y0/v_ph;
Cc = Yc/v_ph_c;
L = 1/v_ph/Y0;
Lc = 1/v_ph_c/Yc;
Ic = 0.0002*thickness/0.00000001*W/0.000001; % samuel's formula
% Ic= 1000000000; % if we only want to consider the coupler's non linearity
Icc = 0.0002*thickness/0.00000001*W_c/0.000001; % samuel's formula
if nockit5_fit
% parameters correction from fit. use these to get somthing close to the
% measurement for 2 traces NOCKIT5, but note that we still have to explain the factor of 2 in
% the phase velocity. the other two factors are close to 1, so they are OK.
x = [1.9935    0.9193    0.8418];  

%x = [2,1,1]
%x =  [1.0479    0.5550    0.4589  -52.1305]
v_ph = v_ph*x(1);
    v_ph_c = v_ph_c*x(2);
    Yc =  Yc/x(3);
    
end

% frequency etc.
freq = [6e9]; 
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
G.Edges.BCval(G.findedge(nodes(input_idx(2),1),nodes(input_idx(2),2))) = pump_amp;


% L_arr = G.Edges.L;
% C_arr = G.Edges.C;
% Len_arr = G.Edges.Len;
% BC_arr = G.Edges.BC;
% BCval_arr = G.Edges.BCval;
% %

%% solve

% which process_graph
% which solve_graph
graph_data = process_graph(G);

[t_edges, r_edges] = solve_graph_non_lin_2(graph_data,freq,iterations(pwr_idx));
%[t_edges, r_edges] = solve_graph(graph_data,freq);

% read solution: this part is specific to the NOCKIT geometry 
 t = reshape(t_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
 r = reshape(r_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);


        trans(pwr_idx,ii,:) = t(:,end);
    end
        trans_norm(pwr_idx,:,:) = trans(pwr_idx,:,:)./trans(pwr_idx,1,:);

 

%plot 
figure(204)
% clf
plot(phase,(abs(trans_norm(pwr_idx,:,2))), 'linewidth', 1.5);
grid on;
xlabel( "phase difference" , "fontsize", 15)
ylabel( "transmission" , "fontsize", 15)
title(sprintf("interference at %g GHz, signal power %g dBm", freq*1e-9, sig_pwr),"fontsize", 15)
% leg = legend(num2str((1:M)'),"location", "best", "fontsize", 16);
%     title(leg, "line")
 set(gca,'XTick',0:pi/2:2*pi) 
 set(gca,'XTickLabel',{'0','\pi/2','\pi','3\pi/2','2\pi'})
 
 
hold on 
end

%%
 cmap = colormap(jet(N_pwrs)) ; %Create Colormap
 cbh = colorbar ; %Create Colorbar
 cbh.Ticks = linspace(0, 1,5) ; %Create 8 ticks from zero to 1
 cbh.TickLabels = num2cell(linspace(pump_min, pump_max,5)) ;    %Replace the labels of these 8 ticks with the numbers 1 to 8
cbh.Label.String = "pump power (dBm)";
cbh.Label.FontSize = 16;

