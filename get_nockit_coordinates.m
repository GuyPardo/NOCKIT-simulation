function [coordinates] = get_nockit_coordinates(G,nockit_params)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
% define coordinates for plotting the graph: (this has no effect on the solution)
N = nockit_params.N;
M = nockit_params.M;
nodes_num = G.numnodes;
L0 = nockit_params.L0;
d = nockit_params.d;

x = repmat(L0*(0:N+1), M,1);
y = repmat(d*fliplr((0:M-1)), 1,N+2); % the y coordinates are in the flipped to plot from top to bottom

x = reshape(x, 1,nodes_num);
y = reshape(y, 1,nodes_num);

coordinates = [x;y];

%G
end

