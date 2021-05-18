function [cost] = get_cost_6(nockit_params,freq,data_dB, x, ind_vec)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
% tic
% x = factor changing [t,W,Wc,H,gap_c,lam]


 display(x)
 
G = change_params(nockit_params,x);

 N = nockit_params.N;
 M= nockit_params.M;
 
 

 

  if x(end)<0
      db_offset = x(end);
  else
      db_offset = -51;
  end
 

nodes = reshape(1:M*(N+2),M,N+2 );  


graph_data = process_graph(G);


trans = nan(M,length(freq));
% loop on frequencies
for i=1:length(freq)
    [t_edges, r_edges] = solve_graph(graph_data,freq(i)); % solve
    
  % read solution: (this part is specific to the NOCKIT geometry)   
    %ref(:,i) = r_edges(G.findedge(nodes(:,1),nodes(:,2)));
   [ref(:,i),trans(:,i)] = read_nockit_solution(nockit_params, G, t_edges, r_edges);
    
end
   trans_dB = 20*log10(abs(trans));
   
% size(trans_dB)
% size(data_dB)
    cost = sum(sum(abs(trans_dB(ind_vec,:) +db_offset - data_dB(ind_vec,:)).^2));
   
%     toc;
end

