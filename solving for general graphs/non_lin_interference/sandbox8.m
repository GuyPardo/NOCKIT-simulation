% sandbox
% here we create a nockit2 grpah and solving with the non linear solver
% 

% close all
clearvars
 nockit5_fit=true;
 coplanar_couplers=false;
%% construct graph
% this part constructs the graph representing the 2 traces ladder network
% (NOCKIT) 

iterations =100;
% intereference experiment setup
N_pwrs = 10;
sig_pwr = -55;% -55; % dbm
pump_pwr = -43; % dbm

phase = pi/2;
sig_amp = sqrt(50*10^((sig_pwr/10) - 3))*exp(1i*phase); % in Volts, assuming here it's about 50 ohm.\
pump_amp = sqrt(50*10.^((pump_pwr/10) - 3));


input_idx = [2,1]; % [signal, pump]

tic

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
% Ic=1e9;
Ic = 0.0002*thickness/0.00000001*W_c/0.000001; % samuel's formula

Ic = 60e-6*(thickness/6e-9)*(W/2.3e-6); % % from mikita measurement
Icc = 60e-6*(thickness/6e-9)*(W_c/2.3e-6); % % from mikita measurement


gnd_conductance = 0;
gnd_conductance_c = gnd_conductance*W_c/W;
% 
% C = C*(1+loss*1i);
% Cc = Cc*(1+loss*1i);

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
freq = [5.1e9]; 
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
G.Edges.gnd_cond(G.Edges.Weight==2) = gnd_conductance*ones(sum(G.Edges.Weight==2),1);
G.Edges.gnd_cond(G.Edges.Weight==1) = gnd_conductance_c*ones(sum(G.Edges.Weight==1),1);
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
try
G.Edges.BCval(G.findedge(nodes(input_idx(2),1),nodes(input_idx(2),2))) = pump_amp;
end
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

[t_edges, r_edges] = solve_graph_non_lin_2(graph_data,freq,iterations,true);
%[t_edges, r_edges] = solve_graph(graph_data,freq);

% read solution: this part is specific to the NOCKIT geometry 
 t = reshape(t_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
 r = reshape(r_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);


%trans(:,ii) = t(:,end);


%%
 %% reconstruct physical quantities
% defining coordinates along the lines:
Npoints = 100; %points per segment

x = linspace(0,(N+1)*L0,(N+1)*Npoints);

% pre-allocating
V = zeros(length(x),M); % voltage
I = zeros(length(x),M); % current


% calculating
for j=1:M
    for n=1:N+1
       V(Npoints*(n-1)+1:Npoints*n,j) = t(j,n)*exp(1i*k0*x(1:Npoints)) + r(j,n)*exp(-1i*k0*x(1:Npoints));
       
       I(Npoints*(n-1)+1:Npoints*n,j) = Y0*(t(j,n)*exp(1i*k0*x(1:Npoints)) - r(j,n)*exp(-1i*k0*x(1:Npoints)));
       
    end
end

% calculate power
P = 0.5*real(V.*conj(I));


%%
%correct 

%% plot - colormap

% figure(203)
% clf
% imagesc(transpose(real(P)), "XData",x )
% shading flat
% colorbar
% yticks(1:M)
% 
% title(sprintf("power propagation at %g GHz", freq*1e-9))
% ylabel( "line" , "fontsize", 15)
% xlabel( "position along line (m)" , "fontsize", 15)

colormap jet

%% plot - graphs
figure(204)
clf
plot(x,real(P), 'linewidth', 1.5);
leg = legend(num2str((1:M)'),"location", "best", "fontsize", 13);
grid on;
xlabel( "position along line (m)" , "fontsize", 15)
ylabel( "power (a.u)" , "fontsize", 15)
title(sprintf("power propagation at %g GHz", freq*1e-9),"fontsize", 15)
