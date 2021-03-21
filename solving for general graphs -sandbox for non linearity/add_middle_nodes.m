function [H] = add_middle_nodes(G,edges,N)
%adds intermediate nodes to a graph
%   this function recieves a gaprh object G, for which the nodes have two
%   attributes: X and Y (representing the node coordinates), and the edges
%   can have any number of properties.
%   the input variable edges is a vector of integers representing
%   edge-indices on the graph.
%   the function removes returns a new graph H, which is the same as G but
%   with the following changes:
%   1. all the edges specifyed by the edges variable are removed
%   2. each such edge is replaced by N+1 new edges with N intermediate
%   nodes
%   the new edges have the same attributes as he edge they replace
%   the Node coordinates are distributed equally between the end nodes of
%   the original edge

H=G;
prop_names = G.Edges.Properties.VariableNames;
edge_vars = G.Edges.Variables;
% edgetable = array2table(G.Edges.Variables, 'variablenames',prop_names(2:end))
for i=1:length(edges)
    nodenum = H.numnodes; % num of nodes before we ass the new ones
    endnodes = G.Edges.EndNodes(edges(i),:); % end nodes for the current (ith) edge
    
    % coordinats for new nodes
    x = linspace(G.Nodes.X(endnodes(1)),G.Nodes.X(endnodes(2)), N + 2);
    x = x(2:end-1);
    y = linspace(G.Nodes.Y(endnodes(1)),G.Nodes.Y(endnodes(2)), N + 2);
    y = y(2:end-1);

    % add new nodes
    node_table = table(x',y', 'variablenames', ["X", "Y"]);
    H=H.addnode(node_table);

    edgetable=array2table(repmat(edge_vars(edges(i),:),N+1,1),'variablenames',prop_names(2:end) );

    % add new edges
    H = H.addedge([endnodes(1) , nodenum + (1:N)] , [nodenum + (1:N), endnodes(2)], edgetable);
    
    % remove edge
    H = H.rmedge(endnodes(1), endnodes(2));

end

end

