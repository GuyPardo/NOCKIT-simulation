
function [G] = change_params(nockit_params,x)
%UNTITLED9 Summary of this function goes here



nockit_params_new.N = nockit_params.N; % number of couplers. (= number of unit cells minus 1) 
nockit_params_new.M = nockit_params.M; % number of lines
nockit_params_new.L0 = nockit_params.L0; % length of each line segment
nockit_params_new.d = nockit_params.d; % length of each coupler segment
nockit_params_new.t = nockit_params.t*x(1); % thickness of metal
nockit_params_new.W = nockit_params.W*x(2) ; % width of main traces
nockit_params_new.W_c = nockit_params.W_c*x(3); % width of couplers
nockit_params_new.H = nockit_params.H*x(4); % thickness of dielectric
nockit_params_new.gap_c = nockit_params.gap_c*x(5); % for coplanar couplers, width of the gap between trace and ground.
nockit_params_new.input_idx =nockit_params.input_idx;

G = get_nockit_graph(nockit_params_new);

 lam = x(6); %inductance factor
 
 G.Edges.v_ph(G.Edges.Weight==2) = G.Edges.v_ph(G.Edges.Weight==2)/sqrt(lam);
G.Edges.v_ph(G.Edges.Weight==1) = G.Edges.v_ph(G.Edges.Weight==1)/sqrt(lam);
G.Edges.Y(G.Edges.Weight==1) = G.Edges.Y(G.Edges.Weight==1)/sqrt(lam);
G.Edges.Y(G.Edges.Weight==2) = G.Edges.Y(G.Edges.Weight==2)/sqrt(lam);


end

