function [trans_dB] =get_trans(G,freq, x, dB_offset, indices)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% tic
% display(x)
if length(x)>4

h = x(1); t = x(2); lam = x(3); W = x(4); Wc = x(5);
else
    h = x(1); t = x(2); lam = x(3); W = 1; Wc = 1;
end
    
M = 7;
N = 31;
nodes = reshape(1:M*(N+2),M,N+2 );  
H = G;
H.Edges.v_ph(H.Edges.Weight==2) = H.Edges.v_ph(H.Edges.Weight==2)*sqrt(t*h/lam);
H.Edges.v_ph(H.Edges.Weight==1) = H.Edges.v_ph(H.Edges.Weight==1)*sqrt(t*h/lam);
H.Edges.Y(H.Edges.Weight==1) = H.Edges.Y(H.Edges.Weight==1)/sqrt(lam*h/t/Wc^2);
H.Edges.Y(H.Edges.Weight==2) = H.Edges.Y(H.Edges.Weight==2)/sqrt(lam*h/t/W^2);

graph_data = process_graph(H);


trans = nan(M,length(freq));
% loop on frequencies
for i=1:length(freq)
    [t_edges, ~] = solve_graph(graph_data,freq(i)); % solve
    
  % read solution: (this part is specific to the NOCKIT geometry)   
    %ref(:,i) = r_edges(H.findedge(nodes(:,1),nodes(:,2)));
    trans(:,i) =t_edges(H.findedge(nodes(:,end-1),nodes(:,end))); 
    
end
   trans_dB = 20*log10(abs(trans(indices,:))) +dB_offset;
   

    
   
%     toc;
end