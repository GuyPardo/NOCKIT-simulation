clearvars
addpath('')

addpath(genpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation'))
%

% read_nockit6_data
load NOCKIT5_2traces_data.mat

%%
 % down sampling
 
 % smooth
  smooth_num=1;
  data_smooth_mag = zeros(length(freq),2);
for i=1:2
     data_smooth_mag(:,i) = smooth((abs(data(:,i))),smooth_num );
end
 data_smooth_dB = 20*log10(data_smooth_mag);
 res = 401;
% res = length(freq);
%data_dB_red = nan(2,res);
freq_red = downsample(freq,res);
data_dB_red = transpose(downsample(data_smooth_dB, res));

figure(401); clf;
plot(freq_red, data_dB_red)


%% fitting

  ind_vec = 1:2;

% geometry: and network structure
get_nockit2_params;

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


%%X = [t,W,Wc,H,lam]
% x0 = [ 0.9575   1.0047  0.9312 1.1149 ,1, 0.8605       ];
x0 = [1,1,1,1,.23];
x0 = [ 0.9385    1.0180    1.0302    1.0933    0.2274];
x0 =   [ 0.9385    1.0180    1.40302    2.233    0.7274];
x0 =   [0.8341    0.9311    1.2711    2.2179    0.6633];
 x0 = [0.8341    0.9311    1.2711     1.0139    0.6633];
% x0=[1.0709    1.0232    0.9961    1.0622    0.9548];
  lb = [.2, .2, .2,0.2,0.2];
  ub = [3, 3,3,3,3];
    
A = [];
b = [];
Aeq = [];
beq = [];

%% minimize
fun = @(x) get_cost_22(nockit_params,freq_red, data_dB_red, x, ind_vec);
fun(x0)
%%
% [X,costval] = fmincon(fun,x0, A,b,Aeq,beq,lb,ub,[], options)
 [X,costval] = fminsearchbnd(fun,x0,lb,ub ,options)

%% plot
%X = [t,W,Wc,H,lam]
% % X = [ 0.9575   1.0047  0.9312 1.1149 ,1, 0.8605       ];
% X =  [1.0709    1.0232    0.9961    1.0622    0.9548];
% X = [1,1,1,1,1]
% X = [0.5044    1.0437    1.3803    1.6273    0.8848];
% X = [ 0.9385    1.0180    1.0302    1.0933    0.2274]
% X = [    0.6695    0.2202    0.2368    1.2618    0.2003]
%    X = [ 0.8074    0.6906    0.7212    1.1616    0.2152]
%  X = [ 0.9385    1.0180    1.40302    2.233    0.7274]
%  X = [0.4850    0.2201    0.3139    2.3889    0.4371]
%    X = [0.6316    0.7151    0.9981    2.2361    0.5197]
%    X = [0.8341    0.9311    1.2711    2.2179    0.6633] % very nice
% X = [0.3661    1.5997    2.1244    1.8366    0.2278]
% X = [  0.8670    0.8932    1.2340    1.0613    0.7141];
X  =  [ 0.9299    1.1904    1.5986    1.6693    1.0784];
N = nockit_params.N;
M = nockit_params.M;
G = change_params(nockit_params,X);
% offset = X(end);
graph_data = process_graph(G);


trans = nan(M,length(freq_red));
% loop on frequencies
for i=1:length(freq_red)
    [t_edges, r_edges] = solve_graph(graph_data,freq_red(i)); % solve
    
  
    
    [trans(:,i),ref(:,i)] = read_nockit_solution(nockit_params, G, t_edges, r_edges);
    
end
   trans_dB = 20*log10(abs(trans)) - 21;
%
figure(333); clf;

% ylim('auto')
% legend("4-->1","4-->2","4-->3","4-->4");
 grid on; xlabel('freq (Hz)'); ylabel('dB'); title('transmission')
 hold on
  cc = colororder();

  for i=ind_vec
   
    plot(freq_red, data_dB_red(i,:), 'color',cc(i,:), 'linewidth', 2)
  end

for i=ind_vec
 plot(freq_red, (trans_dB(i,:)), 'linewidth',2 ,'linestyle', '--','color',cc(i,:));
end

