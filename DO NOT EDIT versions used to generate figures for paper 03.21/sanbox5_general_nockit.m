%% 
% written by Guy  2020_10_20
% solves for an arbitrary graph.
% you  have to define the graph as a digraph (directed graph) object.
% each Edge should have the following properties:
% k - wavenumber for this edge
% Y - characteristic admittance for this edge
% L - length of the edge
% BC - boundary conditions specification according to:
                % 0 - do nothing
                % 1 - set t to 0
                % 2 - set r to 0
                % 3 - set t to 1
                % 4 - set r to 1
% if you don't know how to add properties to graph edges see:
% https://www.mathworks.com/help/matlab/math/add-graph-node-names-edge-weights-and-other-attributes.html#AddGraphNodeNamesEdgeWeightsAndOtherAttributesExample-4
clearvars
 nockit5_fit=true;
 coplanar_couplers=true;
%% construct graph
% this part constructs the graph representing the 2 traces ladder network
% (NOCKIT) 



tic


% geometry: and network structure
N=31; % number of couplers. (= number of unit cells minus 1) 
M = 7; % number of lines
L0 = 100e-6; % length of each line segment
d = 27e-6; % length of each coupler segment
t = 8.5e-9;%10e-9;
W = 3e-6;
W_c = 200e-9;
H = 29e-9;%16e-9;
gap_c = 8e-6;



% % geometry: and network structure
% N=30; % number of couplers. (= number of unit cells minus 1) 
% M = 7; % number of lines
% L0 = 100e-6; % length of each line segment
% d = 20e-6; % length of each coupler segment
% t = 8e-9;
% W = 2.3e-6;
% W_c = 300e-9;
% H = 16e-9;
% gap_c = 1.4e-6;
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
%x =  [1.0479    0.5550    0.4589  -52.1305]
v_ph = v_ph*x(1);
    v_ph_c = v_ph_c*x(2);
    Yc =  Yc/x(3);
    
end

 input_idx = [4];

% % geometry: and network structure
% N=30; % number of couplers. (= number of unit cells minus 1) 
% M = 7; % number of lines
% L0 = 100e-6; % length of each line segment
% d = 20e-6; % length of each coupler segment
% input_idx = [4];   % can be more than one.
% % physical parameters: (see the NOCKIT simulation for the way these values were calculated  )
% v_ph = 1.361104539023962e+06; % phase velocity for lines
% v_ph_c =  1.408763793738406e+06; % phase velocity for couplers
% Y0 = 0.020259488114653; % admittance for lines
% Yc = 0.002735070881675; % admittance for couplers



% % geometry: and network structure
% N=31; % number of couplers. (= number of unit cells minus 1) 
% M = 7; % number of lines
% L0 = 100e-6; % length of each line segment
% d = 27e-6; % length of each coupler segment
% t = 8.5e-9;%10e-9;
% W = 3e-6;
% W_c = 200e-9;
% H = 29e-9;%16e-9;
% gap_c = 1.4e-6;



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



L_arr = G.Edges.L;
v_ph_arr = G.Edges.v_ph;
Y_arr = G.Edges.Y;
BC_arr = G.Edges.BC;

%% graph pre-proccesing:

% define cell arrays of edge-ids:
outedges_cell = cell(1,nodes_num);
inedges_cell = cell(1,nodes_num); 
edges_cell = cell(1,nodes_num); 
for i=1:nodes_num
   outedges_cell{i}  = G.outedges(i);
   inedges_cell{i}  = G.inedges(i);
   edges_cell{i} = [ inedges_cell{i}; outedges_cell{i}];
end

%  index translation:
tIdx =  2*(1:edge_num)-1;
rIdx =  2*(1:edge_num);

disp('graph construction and pre-proccessing:')
toc




%% encode equations
tic


% pre allocation:
mat = zeros(2*edge_num);
vec = zeros(2*edge_num,1);
eqn_count = 0; % the number of equation we have encoded



% loop on nodes
for i=1:nodes_num
    outedges = outedges_cell{i}; % edges going from  node i
    inedges = inedges_cell{i}; % edges going out of  node i
    edges = [ inedges; outedges]; % all edges connected to node i

    
    % edges parameters
    % by converntion we only add exp(+-ikL) for the inedges. so we set
    % L_out and k_out to zero.
    L_in = L_arr(inedges);
    L_out = zeros(size(outedges));
    
    k_out = zeros(size(outedges));
    Y_out  = Y_arr(outedges);
    
    k_in = omega*(v_ph_arr(inedges)).^-1;
    Y_in = Y_arr(inedges);
    
    L = [L_in; L_out];
    k = [k_in; k_out];
    Y = [-Y_in; +Y_out];  % the minus sign is to differentiate btw current into the node and out of the node
    
    
    % encode voltage continuity equation for the  ith node: loop on node edges
    for j = 1:length(edges)-1
       edge = edges(j);
       next_edge = edges(j+1);
        % voltage eqn:
        mat(eqn_count+1,tIdx(edge)) = exp(1i*L(j)*k(j)); 
        mat(eqn_count+1,rIdx(edge)) = exp(-1i*L(j)*k(j));
        mat(eqn_count+1,tIdx(next_edge)) = -exp(1i*L(j+1)*k(j+1));
        mat(eqn_count+1,rIdx(next_edge)) = -exp(-1i*L(j+1)*k(j+1));
        eqn_count = eqn_count+1; 

    end
    
    % encode the ith node current equation:
    if length(edges)>1 % current equations are only relevant for nodes connected to at least 2 edges. 
        mat(eqn_count+1, tIdx(edges)) = Y.*exp(1i*k.*L);
        mat(eqn_count+1, rIdx(edges)) = -Y.*exp(-1i*k.*L);
        eqn_count = eqn_count+1;
    end
 
end

% boundary conditions eqns:
indices = find(BC_arr);% find indices where BC_arr is not zero
for i = 1:length(indices)
    ii = indices(i);
    % add the correct eqns based on the conventions:
    % 1 - set t to 0
    % 2 - set r to 0
    % 3 - set t to 1
    % 4 - set r t0 1
    switch BC_arr(indices(i))
        case 1  
            mat(eqn_count+1, tIdx(ii))=1;
            eqn_count = eqn_count+1;
        case 2
            mat(eqn_count+1, rIdx(ii))=1;
            eqn_count = eqn_count+1;
        case 3
            mat(eqn_count+1, tIdx(ii))=1;
            vec(eqn_count+1)=1;
            eqn_count = eqn_count+1;
        case 4
            mat(eqn_count+1, rIdx(ii))=1;
            vec(eqn_count+1)=1;
            eqn_count = eqn_count+1;
    end           
end        

disp('encode eqns:')
toc
%% solve
tic
solution = linsolve(mat,vec);

t_edges = solution(1:2:end-1);
r_edges = solution(2:2:end);

 
 % check energy conservation:
end_nodes_out = logical(cellfun(@(x) length(x), outedges_cell).*~cellfun(@(x) length(x), inedges_cell));
end_nodes_in = logical(~cellfun(@(x) length(x), outedges_cell).*cellfun(@(x) length(x), inedges_cell));
 
end_edges_in = cell2mat(inedges_cell(end_nodes_in));
end_edges_out = cell2mat(outedges_cell(end_nodes_out));

     
sum_in =  sum(abs(t_edges(end_edges_out)).^2) + sum(abs(r_edges(end_edges_in)).^2); 
sum_out =  sum(abs(r_edges(end_edges_out)).^2) + sum(abs(t_edges(end_edges_in)).^2); 
 cond =abs(sum_in-sum_out);
 if cond>1e-14
     warning('energy conservation condition does not hold')
     cond
 end
%

% read solution: this part is specific to the NOCKIT geometry 
 t = reshape(t_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
 r = reshape(r_edges(G.findedge(nodes(:,1:N+1) ,nodes(:,2:N+2) )), M,N+1);
% 
 disp('solve:')
toc
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


%% plot - colormap

figure(203)
clf
imagesc(transpose(real(P)), "XData",x )
shading flat
colorbar
yticks(1:M)

title(sprintf("power propagation at %g GHz", freq*1e-9))
ylabel( "line" , "fontsize", 15)
xlabel( "position along line (m)" , "fontsize", 15)

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