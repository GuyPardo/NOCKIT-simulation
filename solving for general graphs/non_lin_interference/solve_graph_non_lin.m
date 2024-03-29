function [t_edges, r_edges] = solve_graph_non_lin(graph_data,freq)

% read input
    nodes_num = graph_data.node_num;
    edge_num= graph_data.edge_num;
    Len_arr = graph_data.Len_arr;
    L_arr =  graph_data.L_arr;
    C_arr = graph_data.C_arr;
    BC_arr = graph_data.BC_arr; 
    BCval_arr = graph_data.BCval_arr;
    Ic_arr = graph_data.Ic_arr;


    v_ph_arr = 1./sqrt(L_arr.*C_arr);
    Y_arr = sqrt(C_arr./L_arr);

    
   
% solve linearly:
 [t,r] = solve_graph(graph_data, freq);
 
 % calculate current
 I = (t-r).*Y_arr;
 
 % correct L:
 new_graph_data = graph_data;
 new_graph_data.L_arr = graph_data.L_arr.*(1+(I./graph_data.Ic_arr).^2);

 % solve again:
 
 [t_edges,r_edges] = solve_graph(new_graph_data, freq);

end

