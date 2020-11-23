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
% geometry: and network structure
N=5; %number of colums (ndoes)
M = 5; % number of rows (nodes)
L0 = 100e-6; % length of each line segment
d = 20e-6; % length of each coupler segment
input_idx = [3];   % can be more than one.
% physical parameters: (see the NOCKIT simulation for the way these values were calculated  )
v_ph = 1.361104539023962e+06; % phase velocity for lines
v_ph_c =  1.408763793738406e+06; % phase velocity for couplers
Y0 = 0.020259488114653; % admittance for lines
Yc = 0.002735070881675; % admittance for couplers

% frequency etc.
freq = [6e9]; 
omega= 2*pi*freq;
k0 = 2*pi*freq/v_ph; % wavenumber for lines
kc = 2*pi*freq/v_ph_c; % wavenumber for couplers

% define grpah: define an array of nodes with M rows and Ncolumns. the
% nodes are numbered such that nodes 1:M are the first column, M+1:2*M are
% the socind column  etc.
nodes = reshape(1:M*N,M,N );


nodes_x = (0:N-1)*(L0+2*d);
nodes_y = (0:M-1)*(L0+2*d);

[XX,YY] = meshgrid(nodes_x, nodes_y);
nodes_XX=reshape(XX, M*N,1);
nodes_YY=reshape(YY, M*N,1);


G = digraph();
% define horizontal lines
% connect each node (m,n) to (m, n+1), with weight  2.
% by convention weight 2 means a "main" edge and not a coupler.
s = nodes(:,1:N-1);
t = nodes(:,2:N);
w = 2*ones(size(s));
G = G.addedge(s,t,w);

% define vertical lines
% connect each node (m,n) to (n+1,n) with weight 1.

s = nodes(1:M-1, 1:N);
t = nodes(2:M, 1:N);
w = 2*ones(size(s));
G = G.addedge(s,t,w);

G.Nodes.X = nodes_XX;
G.Nodes.Y = nodes_YY;



G.plot('Xdata', nodes_XX, 'Ydata', nodes_YY)
%%

% % G.plot('xdata', x, 'ydata',y, 'linewidth', LWidths);
endnodes = G.Edges.EndNodes;
for i=1:G.numedges
    G = add_middle_nodes(G,endnodes(i,1),endnodes(i,2), 2, [1,2,1]);
end


G.plot('Xdata',G.Nodes.X, 'Ydata', G.Nodes.Y)
% define boundary conditions attribute for each edge according to the
% following convention:
% 0 - do nothing
% 1 - set t to 0
% 2 - set r to 0
% 3 - set t to 1
% 4 - set r to 1
%%
G.Edges.BC = zeros(G.numedges,1);
 G.Edges.ID = (1:G.numedges)';
plot(G, 'edgelabel', G.Edges.ID, 'xdata', G.Nodes.X, 'Ydata', G.Nodes.Y)

% add inputs/outputs
% all inputs and outputs are represented by OUT-going edges.


% left :  
% left BC : all r's are zero exctept the special one(s): (it's an OUT going edge, so r is going to the right)
BC_left_arr = 2*ones(M,1);
BC_left_arr(input_idx) = 4;
edge_table  = table(ones(M,1),BC_left_arr , G.numedges + (1:M)','variablenames', {'Weight', 'BC','ID'});
[G, left_edges] = add_free_edge(G,nodes(:,1), edge_table, -L0,0);
left_edges_ID = G.Edges.ID(left_edges);
plot(G, 'edgelabel',G.Edges.ID);

% right: 
 edge_table  = table(ones(M,1),2*ones(M,1),G.numedges + (1:M)' ,'variablenames', {'Weight', 'BC','ID'});
 [G, right_edges] = add_free_edge(G,nodes(:,N), edge_table, L0,0);
 right_edges_ID = G.Edges.ID(right_edges);

plot(G, 'edgelabel', 1:G.numedges);

% top: 
edge_table  = table(ones(N,1),2*ones(N,1),G.numedges + (1:N)' ,'variablenames', {'Weight', 'BC', 'ID'});
[G, top_edges] = add_free_edge(G,(nodes(1,:))', edge_table,0,-L0);
plot(G, 'edgelabel', G.Edges.ID)
top_edges_ID = G.Edges.ID(top_edges);


% bottom: 
edge_table  = table(ones(N,1),2*ones(N,1), G.numedges + (1:N)','variablenames', {'Weight', 'BC', 'ID'});
[G, bottom_edges]  = add_free_edge(G,nodes(M,:), edge_table,0,L0);
bottom_edges_ID = G.Edges.ID(bottom_edges);


edgeOrder(G.Edges.ID) = 1:G.numedges;

%
% define edges attributres; phase velocity, length and characteristic
% admittance:
% rememeber weight = 1 means coupler edge, weight = 2 means regular edge:
% clearvars G.Edges.v_ph G.Edges.L G.Edges.Y
G.Edges.v_ph(G.Edges.Weight==2) = v_ph*ones(sum(G.Edges.Weight==2),1);
G.Edges.v_ph(G.Edges.Weight==1) = v_ph_c*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.Y(G.Edges.Weight==2) = Y0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Y(G.Edges.Weight==1) = Yc*ones(sum(G.Edges.Weight==1),1);






L_arr = G.Edges.L;
v_ph_arr = G.Edges.v_ph;
Y_arr = G.Edges.Y;
BC_arr = G.Edges.BC;

clf
G.plot('edgelabel',G.Edges.ID,'layout', 'force')
%% graph pre-proccesing:
edge_num = G.numedges;
nodes_num = G.numnodes;
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

 
 % check energy conservation: (for a general graph!)
end_nodes_out = logical(cellfun(@(x) length(x), outedges_cell).*~cellfun(@(x) length(x), inedges_cell));
end_nodes_in = logical(~cellfun(@(x) length(x), outedges_cell).*cellfun(@(x) length(x), inedges_cell));
 
end_edges_in = cell2mat(inedges_cell(end_nodes_in));
end_edges_out = cell2mat(outedges_cell(end_nodes_out));

     
sum_in =  sum(abs(t_edges(end_edges_out)).^2) + sum(abs(r_edges(end_edges_in)).^2); 
sum_out =  sum(abs(r_edges(end_edges_out)).^2) + sum(abs(t_edges(end_edges_in)).^2); 
 cond =abs(sum_in-sum_out);
 if cond>1e-16
     warning('energy conservation condition does not hold')
     cond
 end
%


% 
 disp('solve:')
toc
%%
% define edge power value
G.Edges.power = abs(t_edges).^2 - abs(r_edges).^2;
 %% plot
 reverse_idx = G.Edges.power<0;
 ids = 1:G.numedges;
 G_rev = G.flipedge(ids(reverse_idx));
 
h = plot(G_rev,'linewidth', 8,  'arrowsize', 25 ,'edgealpha',1 ,'xdata',G_rev.Nodes.X ,'ydata',G_rev.Nodes.Y);
colormap default;
h.EdgeCData = abs(G_rev.Edges.power);
colorbar
 
 