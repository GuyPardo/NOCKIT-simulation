function [trans,ref] = read_nockit_solution(nockit_params, G, t_edges, r_edges)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here
N = nockit_params.N;
M = nockit_params.M;
nodes = reshape(1:M*(N+2),M,N+2 );  

    ref = r_edges(G.findedge(nodes(:,1),nodes(:,2)));
    trans =t_edges(G.findedge(nodes(:,end-1),nodes(:,end)));


end

