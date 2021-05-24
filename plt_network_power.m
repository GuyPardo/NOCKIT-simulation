function [] = plt_network_power(G,coordinates,t_edges,r_edges)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here

Y = sqrt(G.Edges.C./G.Edges.L);

G.Edges.power = Y.*(abs(t_edges).^2 - abs(r_edges).^2);

 reverse_idx = G.Edges.power<0;
 ids = 1:G.numedges;
 G_rev = G.flipedge(ids(reverse_idx));
 figure(gcf);
h = plot(G_rev,'linewidth', 4,  'arrowsize', 12 ,'edgealpha',1 ,'xdata',coordinates(1,:) ,'ydata',coordinates(2,:));
colormap default;
h.EdgeCData = abs(G_rev.Edges.power);
colorbar 
%G.plot('xdata', coordinates(1,:), 'ydata',coordinates(2,:), 'linewidth', LWidths);
end

