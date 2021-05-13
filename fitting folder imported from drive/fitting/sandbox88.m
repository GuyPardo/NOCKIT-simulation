clearvars
addpath('')

addpath(genpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation'))
%%
read_nockit6_data
%%
 % down sampling
res = 201;
%data_dB_red = nan(2,res);
freq_red = downsample(freq,res);
data_dB_red = (downsample(data_smooth_dB', res))';


%% fitting



% geometry: and network structure
nockit_params.N=31; % number of couplers. (= number of unit cells minus 1) 
nockit_params.M = 7; % number of lines
nockit_params.L0 = 100e-6; % length of each line segment
nockit_params.d = 27e-6; % length of each coupler segment

nockit_params.W=3e-6; % width of primary and secondary transmission lines
nockit_params.t=8.5e-9; % thickness of WSi (sputtered)
nockit_params.H=35e-9; % height of dielectric (say, Si - evaporated)
nockit_params.W_c=200e-9; % width of coupling line
nockit_params.gap_c = 8e-6;
nockit_params.input_idx = [4];   % can be more than one.

G = get_nockit_graph(nockit_params);


% frequency etc.
%freq = 1e9*linspace(3,9,201); 
omega= 2*pi*freq;





% graph_data = process_graph(G);


%% config
% options = optimoptions('fmincon','Display','iter', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% 
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');

options = optimset('PlotFcns',@optimplotfval, 'Display','iter');


%%X = [t,W,Wc,H,gap_c,lam]
x0 = [1,1,1,1,1,1]; 

  lb = [.6, .6, .6,0.6,0.8,.2];
  ub = [1.3, 1.3,1.3,1.3,1.2,2];
    
A = [];
b = [];
Aeq = [];
beq = [];

%% minimize
fun = @(x) get_cost_6(nockit_params,freq_red, data_dB_red, x);
fun(x0)
%%
% [X,costval] = fmincon(fun,x0, A,b,Aeq,beq,lb,ub,[], options)
 [X,costval] = fminsearchbnd(fun,x0,lb,ub ,options)

%% plot
G = change_params(nockit_params,X);
% offset = X(end);
graph_data = process_graph(G);


trans = nan(M,length(freq));
% loop on frequencies
for i=1:length(freq)
    [t_edges, ~] = solve_graph(graph_data,freq(i)); % solve
    
  % read solution: (this part is specific to the NOCKIT geometry)   
    %ref(:,i) = r_edges(G.findedge(nodes(:,1),nodes(:,2)));
    trans(:,i) =t_edges(G.findedge(nodes(:,end-1),nodes(:,end))); 
    
end
   trans_dB = 20*log10(abs(trans)) - 51;
%%
figure(333); clf;

% ylim('auto')
% legend("4-->1","4-->2","4-->3","4-->4");
 grid on; xlabel('freq (Hz)'); ylabel('dB'); title('transmission')
 hold on
  cc = colororder();
  
  for i=2:4
   
    plot(freq, data_smooth_dB(i,:), 'color',cc(i,:), 'linewidth', 2)
  end

for i=2:4
 plot(freq, (trans_dB(i,:)), 'linewidth',2 ,'linestyle', '--','color',cc(i,:));
end

