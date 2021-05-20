function [graph_data] = process_graph_NL(G)
% written by guy 2021_04_21. pre processing of a directed graph G
% this function recives a digraph object G and returns a struct graph_data
% with the following fields:
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

% 

    graph_data.node_num = G.numnodes;
    graph_data.edge_num = G.numedges;
    
    graph_data.Len_arr = G.Edges.Len;
    graph_data.L_arr = G.Edges.L;
    graph_data.C_arr = G.Edges.C;
    graph_data.BC_arr = G.Edges.BC;
    graph_data.BCval_arr = G.Edges.BCval;
    graph_data.Ic_arr = G.Edges.Ic;
    graph_data.gnd_cond_arr = G.Edges.gnd_cond;
    % define cell arrays of edge-ids:
    graph_data.outedges_cell = cell(1,graph_data.node_num);
    graph_data.inedges_cell = cell(1,graph_data.node_num); 
    for i=1:graph_data.node_num
       graph_data.outedges_cell{i}  = G.outedges(i);
       graph_data.inedges_cell{i}  = G.inedges(i);
    end
    
   


end