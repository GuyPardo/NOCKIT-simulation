clearvars

nockit5_fit=false;
coplanar_couplers=true;
freq = 6e9;
omega = 2*pi*freq;
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
H = 35e-9;
gap_c = 8e-6;
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


input_idx = [4];   % can be more than one.
% critical current: formulas taken grom Samuel's microstrip calculator
Ic = 0.0002*t/0.00000001*W/0.000001;
Ic_c = 0.0002*t/0.00000001*W_c/0.000001;


L = 1/v_ph/Y0; Lc = 1/v_ph_c/Yc;
C = 1/v_ph*Y0; Cc = 1/v_ph_c*Yc;

%
%

%%
nodes = reshape(1:M*(N+2),M,N+2 );  
G = digraph();
% define main lines:
% connect each node (m,n) to (m, n+1), with weight  2.
% by convention weight 2 means a "main" edge and not a coupler.
source = nodes(:,1:N+1);
target = nodes(:,2:N+2);
w = 2*ones(size(source));
G = G.addedge(source,target,w);

% define coupling lines:
% connect each node (m,n) to (n+1,n) with weight 1.
% by convention weight 1 means a coupling edge.
source = reshape(nodes(1:M-1, 2:N+1), 1,(N)*(M-1));
target = reshape(nodes(2:M, 2:N+1), 1,(N)*(M-1));
w = ones(size(source));
G = G.addedge(source,target,w);


edge_num = G.numedges;
nodes_num = G.numnodes;

% define coordinates for plotting the graph: (this has no effect on the solution)
x = repmat(L0*(0:N+1), M,1);
y = repmat(d*fliplr((0:M-1)), 1,N+2); % the y coordinates are in the flipped to plot from top to bottom

G.Nodes.X = reshape(x, 1,nodes_num)';
G.Nodes.Y = reshape(y, 1,nodes_num)';
LWidths = 6*G.Edges.Weight/max(G.Edges.Weight); 
figure(504)
clf
G.plot('xdata', G.Nodes.X, 'ydata',G.Nodes.Y, 'linewidth', LWidths);


%%


% define edges attributres; inducrtance (L), and capacitalce (C) per unoit length 
% rememeber weight = 1 means coupler edge, weight = 2 means regular edge:
clearvars G.Edges.C G.Edges.L  G.Edges.Ic
G.Edges.L(G.Edges.Weight==2) = L*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = Lc*ones(sum(G.Edges.Weight==1),1);
G.Edges.C(G.Edges.Weight==2) = C*ones(sum(G.Edges.Weight==2),1);
G.Edges.C(G.Edges.Weight==1) = Cc*ones(sum(G.Edges.Weight==1),1);
G.Edges.Ic(G.Edges.Weight==2) = Ic*ones(sum(G.Edges.Weight==2),1);
G.Edges.Ic(G.Edges.Weight==1) = Ic_c*ones(sum(G.Edges.Weight==1),1);


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



%% solve

graph_data = process_graph2(G);
[t_edges, r_edges] = solve_graph_non_lin(graph_data,freq, 1);

%% define new graph
subdivide = 3;
H = add_middle_nodes(G,(1:G.numedges), subdivide - 1);

figure(45345)
clf
H.plot('xdata', H.Nodes.X, 'ydata',H.Nodes.Y);

%% 
for i=1:G.numedges

L = G.L(i);
C = G.C(i);
RR = 0;
GG = 1e-6;


    k  = 1i*sqrt((RR - 1i*omega*L)*(GG - 1i*omega*C));
   Y = sqrt((GG - 1i*omega*C)/(RR - 1i*omega*L));
    
    pos = H.Nodes.
     I = t_edges(i)*exp(1i*k*)
end



%%




% read solution: this part is specific to the NOCKIT geometry 
 t = reshape(t_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
 r = reshape(r_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);







