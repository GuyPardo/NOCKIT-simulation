function [G, edges] = add_free_edge(G,node_id, attributes)
% written by guy 2020_11_05 adds a "free edge" connected only to one node.
% In matlab graph and digraph objects there is no such thing: each edge has
% to have 2 nodes, so what we do in this function is to create a new node
% and a new edge. attributes is optional parameter (however note that if G
% has weights and or attributes it is required) used to define edge
% properties and/or weight. it should be either a vector of wheights, the
% same size as node_id, or a table according to the matlab conventions. see
% for example :
% https://www.mathworks.com/help/matlab/ref/graph.addedge.html#buofhkr-1-EdgeTable
% returns the new modified  graph (G) and an array of the new edges
% indices (edges)


N = numel(node_id);
previous_numnodes =  G.numnodes;
if nargin==3
    G=G.addedge(node_id, (1:N) + previous_numnodes, attributes);
else
    G = G.addedge(node_id, (1:N) + previous_numnodes);
end


edges = G.findedge(node_id, (1:N) + previous_numnodes);

end