function [t_edges, r_edges] = solve_graph_non_lin_2(graph_data,freq, iterations,plot_iterations)
if nargin<4
    plot_iterations=false;
end
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
    k_arr = 2*pi*freq*(v_ph_arr).^-1;

     % solve linearly:
 [t,r] = solve_graph(graph_data, freq);
 
 for i=1:iterations   
   

 t_previous = t;
 r_previous = r;
 % calculate current (average ebggining and end of each segment)
 I  = real((t-r).*Y_arr);
 I_avr = real(((t-r).*Y_arr  + Y_arr.*(t.*exp(1i*k_arr.*Len_arr)  - r.*exp(-1i*k_arr.*Len_arr)))/2);
%   I = real((t-r).*Y_arr);
 
 % correct L:
 graph_data.L_arr = L_arr.*(1+(I_avr./graph_data.Ic_arr).^2);

 %solve again
 [t,r] = solve_graph(graph_data, freq);
 
 %calculate difference:
 
 diff_r = r-r_previous;
 diff_t = t-t_previous;
 
 if plot_iterations
     figure(54)
    colororder(jet(iterations))
     plot(I)
     hold on
 end
 
 diff_all(i) = sum(diff_r) + sum(diff_t);
 sum_all(i) = sum(r) + sum(t); 
 end
 

 if plot_iterations
figure(568)
subplot(2,1,1);  plot(1:iterations, abs(diff_all)); title("diffrences"); xlabel("iterations");
 subplot(2,1,2); plot(1:iterations, abs(sum_all)); title("values"); xlabel("iterations");
 
 end
 t_edges = t;
 r_edges = r;
 
end

