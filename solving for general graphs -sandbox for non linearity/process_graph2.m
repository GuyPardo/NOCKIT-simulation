function [graph_data] = process_graph2(G)
% written by guy 2020_10_27. pre processing of a directed graph G
% this function recives a digraph object G and returns a struct graph_data
% with the following fields:
    %     nodes_num = number of nodes in the graph edge_num= gnumber of edges
    %     on the graph L_arr = a vector of length edge_num storing the lengths
    %     of the edges v_ph_arr =  a vector of length edge_num storing the
    %     phase valocities of the edges Y_arr = a vector of length edge_num
    %     storing the characteristic admittances of the edges BC_arr = a vector
    %     of length edge_num storing the boundary condition value of each edges
    %     according to the following convetion:
            % 1 - set t to 0
            % 2 - set r to 0
            % 3 - set t to 1
            % 4 - set r t0 1
    %    outedges_cell = a cell array of length node_num where the ith cell is a
    %    collum vector containing the indices of the out-going edges of the ith
    %    nodes
    %    inedges_cell = similarly for in going edges.
% 

    graph_data.node_num = G.numnodes;
    graph_data.edge_num = G.numedges;
    
    
    graph_data.v_ph_arr = sqrt(1./(G.Edges.L.*G.Edges.C));
    graph_data.Y_arr = sqrt(G.Edges.C./G.Edges.L);
    graph_data.BC_arr = G.Edges.BC;
    
    % get length of segments from end nodes:
    node1 = G.Edges.EndNodes(:,1);
    node2 = G.Edges.EndNodes(:,2);
    
    graph_data.L_arr = sqrt((G.Nodes.X(node1) - G.Nodes.X(node2)).^2  + (G.Nodes.Y(node1) - G.Nodes.Y(node2)).^2);
    
    % define cell arrays of edge-ids:
    graph_data.outedges_cell = cell(1,graph_data.node_num);
    graph_data.inedges_cell = cell(1,graph_data.node_num); 
    for i=1:graph_data.node_num
       graph_data.outedges_cell{i}  = G.outedges(i);
       graph_data.inedges_cell{i}  = G.inedges(i);
    end
    
   


end