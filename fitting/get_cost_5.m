function [cost] = get_cost_5(G,freq,data_dB, x)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% tic
% x = factor changing [, H,t,lamnda]
 display(x)
 if length(x)>4

H = x(1); t = x(2); lam = x(3); W = x(4); Wc = x(5);
else
    H = x(1); t = x(2); lam = x(3); W = 1; Wc = 1;
end
    
 
M = 2;
N = 30;
nodes = reshape(1:M*(N+2),M,N+2 );  

G.Edges.v_ph(G.Edges.Weight==2) = G.Edges.v_ph(G.Edges.Weight==2)*sqrt(t*H/lam);
G.Edges.v_ph(G.Edges.Weight==1) = G.Edges.v_ph(G.Edges.Weight==1)*sqrt(t*H/lam);
G.Edges.Y(G.Edges.Weight==1) = G.Edges.Y(G.Edges.Weight==1)/sqrt(lam*H/t/Wc^2);
G.Edges.Y(G.Edges.Weight==2) = G.Edges.Y(G.Edges.Weight==2)/sqrt(lam*H/t/W^2);

graph_data = process_graph(G);


trans = nan(M,length(freq));
% loop on frequencies
for i=1:length(freq)
    [t_edges, ~] = solve_graph(graph_data,freq(i)); % solve
    
  % read solution: (this part is specific to the NOCKIT geometry)   
    %ref(:,i) = r_edges(G.findedge(nodes(:,1),nodes(:,2)));
    trans(:,i) =t_edges(G.findedge(nodes(:,end-1),nodes(:,end))); 
    
end
   trans_dB = 20*log10(abs(trans));
   
% size(trans_dB)
% size(data_dB)
    cost = sum(sum(abs(trans_dB(:,:) -22 - data_dB(:,:)).^2));
   
%     toc;
end

