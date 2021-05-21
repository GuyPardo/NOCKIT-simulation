function [t_edges, r_edges] = solve_graph_NL_envelope(graph_data,freq, iterations,average_current,average_result,plot_iterations)
% if nargin<4
%     average_current = 1;
% end
% if nargin<5
%     plot_iterations=false;
% end
% read input
    nodes_num = graph_data.node_num;
    edge_num= graph_data.edge_num;
    Len_arr = graph_data.Len_arr;
    L_arr =  graph_data.L_arr;
    C_arr = graph_data.C_arr;
    BC_arr = graph_data.BC_arr; 
    BCval_arr = graph_data.BCval_arr;
    Ic_arr = graph_data.Ic_arr;
    I_star_arr = Ic_arr;

    v_ph_arr = 1./sqrt(L_arr.*C_arr);
    Y_arr = sqrt(C_arr./L_arr);
    k_arr = 2*pi*freq*(v_ph_arr).^-1;

    
    
     % solve linearly:
 [t,r] = solve_graph_NL(graph_data, freq);
 t_sum = zeros(1,edge_num);
 r_sum = zeros(1,edge_num);
 
 I_iter = zeros(edge_num,iterations);
 for i=1:iterations   
   
    
 t_previous = t;
 r_previous = r;
 % calculate current (average ebggining and end of each segment)
 Iin  = real((t-r).*Y_arr);
 Iout = real(Y_arr.*(t.*exp(1i*k_arr.*Len_arr)  - r.*exp(-1i*k_arr.*Len_arr)));
 I_mid = (Iin+Iout)/2; % average of two ends of edge
%   I = real((t-r).*Y_arr);
 I_iter(:,i) = I_mid;
 % correct L:
 
 if i>=average_current
    I_avr = sum(I_iter(:,(i-average_current+1):i),2)/average_current;
 else
    I_avr = I_mid;
 end
%  disp(" size I_avr")
%   size(I_avr)
 graph_data.L_arr = L_arr.*(1+(I_avr./I_star_arr).^2);

 %solve again
 [t,r] = solve_graph_NL(graph_data, freq);
 
 %average result
 start_average = iterations - average_result +1;
 if i>start_average
     t_sum = t_sum + t;
     r_sum = r_sum + r;
 end
 
 %calculate difference:
 
 diff_r = r-r_previous;
 diff_t = t-t_previous;
 
 if plot_iterations
%      figure(54)
%     colororder(jet(iterations))
%     subplot(2,1,1)
%      plot(Iin)
%      hold on
%      subplot(2,1,2)
%      plot(Iout)
%      hold on
 end
 
 diff_all(i) = sum(diff_r) + sum(diff_t);
 sum_all(i) = sum(r) + sum(t); 
 
     if i>5
     stop_cond  = all(abs(diff_all(i-4:i)) < 3 e-2*max(abs(diff_all(1:i))));
         if stop_cond
             stop_idx = i;
             break
         else
             stop_idx = iterations;

         end
     end
 end
 

 if plot_iterations
figure(568)
subplot(2,1,1);  plot(1:stop_idx, abs(diff_all)); title("diffrences"); xlabel("iterations");
 subplot(2,1,2); plot(1:stop_idx, abs(sum_all)); title("values"); xlabel("iterations");
 
 end
 
 t_edges = t_sum/average_result;
 r_edges = r_sum/average_result;
 
 if average_result ==1
     t_edges = t;
     r_edges = r;
 end
%  t_edges = t;
%  r_edges = r;
%  
end

