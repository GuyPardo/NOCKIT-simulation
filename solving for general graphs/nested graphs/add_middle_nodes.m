function G = add_middle_nodes(G, node1, node2,N, w)
%written by Guy 2020_10_2
% for two connected nodes node1 and node2, this function adds N
% intermediate nodes.
 % more explicitely - what is done is the edge brw node1 and node2 is
 % deleted, N new nodes are created, and then N+1 new  edges are
 % added instead of the deleted adge. 
 % w is an optional parameter indicating the weigths of the new edges 
 % it should be a vector of length N+1
 
nodes_num = G.numnodes; % initial number of nodes
G = G.addnode(N);
G = G.rmedge( node1, node2); % remove the edge btw node1 and node2

source = [node1 , nodes_num+1:nodes_num+N];
target = [ nodes_num+1:nodes_num+N, node2];
if nargin==5
    G = G.addedge(source, target, w);
else
    G = G.addedge(source, target);
end






end