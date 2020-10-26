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

%% construct graph
% this part constructs the graph representing the 2 traces ladder network
% (NOCKIT) 


clearvars
tic
%geometry:
N=30; % number of couplers. (= number of unit cells minus 1) 
L0 = 100e-6; % length of each line segment
d = 20e-6; % length of each coupler segment
% physical parameters: (see the NOCKIT simulation  )
v_ph = 1.361104539023962e+06; % phase velocity for lines
v_ph_c =  1.408763793738406e+06; % phase velocity for couplers
Y0 = 0.020259488114653; % admittance for lines
Yc = 0.002735070881675; % admittance for couplers

freq = 6e9; 
omega= 2*pi*freq;
k0 = 2*pi*freq/v_ph; % wavenumber for lines
kc = 2*pi*freq/v_ph_c; % wavenumber for couplers

% define ladder grpah:
% the first line is odd nodes, the second line is even nodes.
% such that nodes 1 and 2 are "inputs" and nodes 2*N+3 and 2*N+4 are "outputs"
source = [3:2:2*N+1, 1:2:2*N+1, 2:2:2*N+2];
target = [4:2:2*N+2, 3:2:2*N+3, 4:2:2*N+4];
w = [1*ones(1,N), 2*ones(1,length(source)-N)]; % 1 is for coupling edges, 2 for main edges

G = digraph(source,target,w);
edge_num = G.numedges;
nodes_num = G.numnodes;
clearvars G.Edges.k G.Edges.L G.Edges.Y
G.Edges.v_ph(G.Edges.Weight==2) = v_ph*ones(sum(G.Edges.Weight==2),1);
G.Edges.v_ph(G.Edges.Weight==1) = v_ph_c*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.Y(G.Edges.Weight==2) = Y0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Y(G.Edges.Weight==1) = Yc*ones(sum(G.Edges.Weight==1),1);

% boundary conditions:
% 0 - do nothing
% 1 - set t to 0
% 2 - set r to 0
% 3 - set t to 1
% 4 - set r t0 1
G.Edges.BC = zeros(edge_num,1);
G.Edges.BC(G.findedge(1,3)) = 3;
G.Edges.BC(G.findedge(2,4)) = 1;
G.Edges.BC(G.findedge(2*N+1,2*N+3)) = 2;
G.Edges.BC(G.findedge(2*N+2,2*N+4)) = 2;


L_arr = G.Edges.L;
v_ph_arr = G.Edges.v_ph;
Y_arr = G.Edges.Y;
BC_arr = G.Edges.BC;




figure(504)
G.plot('EdgeLabel', G.Edges.Weight);
outedges_cell = cell(1,nodes_num);
inedges_cell = cell(1,nodes_num);
edges_cell = cell(1,nodes_num);

for i=1:nodes_num
   outedges_cell{i}  = G.outedges(i);
   inedges_cell{i}  = G.inedges(i);
   edges_cell{i} = [ inedges_cell{i}; outedges_cell{i}];
end

disp('graph construction:')
toc




%% encode equations
tic



%  index translation functions:
tIdx =  2*(1:edge_num)-1;
rIdx =  2*(1:edge_num);

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
    L_in = L_arr(inedges);
    L_out = zeros(size(outedges));
    
    k_out = omega*(v_ph_arr(outedges)).^-1;
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
for i = 1:edge_num
    switch BC_arr(i)
        case 1  
            mat(eqn_count+1, tIdx(i))=1;
            eqn_count = eqn_count+1;
        case 2
            mat(eqn_count+1, rIdx(i))=1;
            eqn_count = eqn_count+1;
        case 3
            mat(eqn_count+1, tIdx(i))=1;
            vec(eqn_count+1)=1;
            eqn_count = eqn_count+1;
        case 4
            mat(eqn_count+1, rIdx(i))=1;
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

% read solution: this part is specific to the ladder network
 t = [transpose(t_edges(G.findedge(1:2:2*N+1,3:2:2*N+3))); transpose(t_edges(G.findedge(2:2:2*N+2,4:2:2*N+4)))];
 r = [transpose(r_edges(G.findedge(1:2:2*N+1,3:2:2*N+3))); transpose(r_edges(G.findedge(2:2:2*N+2,4:2:2*N+4)))];
 
 % check energy conservation (this is also specific to the ladder network)
 cond =abs(-1+ abs(r(1,1))^2 + abs(r(2,1))^2+ abs(t(1,end))^2 + abs(t(2,end))^2);
 if cond>1e14
     warning('energy conservation condition does not hold')
     cond
 end
disp('solve:')
toc
 %% reconstruct physical quantities
% defining coordinates along the lines:
Npoints = 100; %points per segment
M=2;
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


%% plot graphs
figure(502)
plot(x,P, 'linewidth', 1.5); grid on; legend('line 1', 'line 2','fontsize',14, 'location', 'best');

xlabel( "position along line (m)" , "fontsize", 15)
ylabel( "power (a.u.)" , "fontsize", 15)
title(sprintf("power propagation at %g GHz", freq*1e-9),"fontsize", 15)
%% plot - colormap
figure(503)
clf
imagesc(transpose(real(P)), "XData",x )
shading flat
colorbar
yticks(1:M)

title(sprintf("power propagation at %g GHz", freq*1e-9))
ylabel( "line" , "fontsize", 15)
xlabel( "position along line (m)" , "fontsize", 15)

colormap jet