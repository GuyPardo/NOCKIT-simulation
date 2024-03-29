function [t_edges, r_edges] = solve_graph(graph_data,freq)
% written by Guy 2020_10_27
% solves for a given graph  and a given frequency
% 
% input arguments:
% graph_data: a struct with specific structure that is
%           generated by running : graph_data = process_graph(G) where G is some
%           digraph object alternatively, you can build the struct manually: the
            % struct should have the following fields:
                %     nodes_num = number of nodes in the graph
                %     edge_num= gnumber of edges on the graph
                %     Len_arr = a vector of length edge_num storing the lengths of the edges
                %     L_arr =  a vector of length edge_num storing the inductances edges
                %     C_arr = a vector of length edge_numstoring the capacitances of the edges
                %     BC_arr = a vector of length edge_num storing the boundary condition value of each edges
                %     according to the following convetion:
                        % 1 - set t to 0
                        % 2 - set r to 0
                        % 3 - set t to the corresponiding BCval
                        % 4 - set r t0 the corresponiding BCval
%                        
                %    BCval_arr  : a vector of length edge_num storing the values for the boudary conditions 
                %    outedges_cell = a cell array of length node_num where the ith cell is a
                %    collum vector containing the indices of the out-going edges of the ith
                %    nodes
                %    inedges_cell = similarly for in going edges.
% freq: a scalar. the frequecy for which to solve. 
% output arguments:
% t_edges : a vector of length graph_data.num_edges containing the computed t
%           amplitudes for each edge
%  r_edges : a vector of length graph_data.num_edges containing the  computed r
%            amplitudes for each edge


%% read input
    nodes_num = graph_data.node_num;
    edge_num= graph_data.edge_num;
    Len_arr = graph_data.Len_arr;
    L_arr =  graph_data.L_arr;
    C_arr = graph_data.L_arr;
    BC_arr = graph_data.BC_arr; 
    BCval_arr = graph_data.BCval_arr;
    outedges_cell = graph_data.outedges_cell; 
    inedges_cell = graph_data.inedges_cell; 

    
v_ph_arr = 1./sqrt(L_arr.*C_arr);
Y_arr = sqrt(C_arr./L_arr);
    
   omega = freq*2*pi; % angular frequency
   % index translation
   tIdx =  2*(1:edge_num)-1;
   rIdx =  2*(1:edge_num);
         
%% encode equations


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
    Len_in = Len_arr(inedges);
    Len_out = zeros(size(outedges));
    
    k_out = zeros(size(outedges));
    Y_out  = Y_arr(outedges);
    
    k_in = omega*(v_ph_arr(inedges)).^-1;
    Y_in = Y_arr(inedges);
    
    Len = [Len_in; Len_out];
    k = [k_in; k_out];
    Y = [-Y_in; +Y_out];  % the minus sign is to differentiate btw current into the node and out of the node
    
    
    % encode voltage continuity equation for the  ith node: loop on node edges
    for j = 1:length(edges)-1
       edge = edges(j);
       next_edge = edges(j+1);
        % voltage eqn:
        mat(eqn_count+1,tIdx(edge)) = exp(1i*Len(j)*k(j)); 
        mat(eqn_count+1,rIdx(edge)) = exp(-1i*Len(j)*k(j));
        mat(eqn_count+1,tIdx(next_edge)) = -exp(1i*Len(j+1)*k(j+1));
        mat(eqn_count+1,rIdx(next_edge)) = -exp(-1i*Len(j+1)*k(j+1));
        eqn_count = eqn_count+1; 

    end
    
    % encode the ith node current equation:
    if length(edges)>1 % current equations are only relevant for nodes connected to at least 2 edges. 
        mat(eqn_count+1, tIdx(edges)) = Y.*exp(1i*k.*Len);
        mat(eqn_count+1, rIdx(edges)) = -Y.*exp(-1i*k.*Len);
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
            vec(eqn_count+1)=BCval_arr(ii);
            eqn_count = eqn_count+1;
        case 4
            mat(eqn_count+1, rIdx(ii))=1;
            vec(eqn_count+1)=BCval_arr(ii);
            eqn_count = eqn_count+1;
    end           
end        



%% solve

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
 if cond>1e-13
     warning('energy conservation condition does not hold:')
     cond
 end





end