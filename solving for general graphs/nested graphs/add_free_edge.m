function [G, edges] = add_free_edge(G,node_id, attributes,dx,dy)
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
X = G.Nodes.X(node_id)+dx;
Y = G.Nodes.Y(node_id) + dy;
nodeprops = table(X,Y, 'variablenames', {'X', 'Y'});

G = G.addnode( nodeprops);


G = G.addnode(N);

    G=G.addedge(node_id, (1:N) + previous_numnodes, attributes);



edges = G.findedge(node_id, (1:N) + previous_numnodes);

end