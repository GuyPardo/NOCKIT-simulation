clearvars
addpath('')

addpath(genpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation'))

name = "\Transmission@-40";

lab = false;


if lab

   folder = 'Z:\Measurements\2020.dir\01_12_20 cooldown\Nockit6';
else
   folder = 'C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\nockit6 data\Nockit6'; 
end
%% import data
smoothing = 0;

% traceid = 1:4;
% path =  strcat(folder,name);
% clearvars data
% for i=1:length(traceid)
% 
%         foldername = sprintf("4to%s", num2str(traceid(i)));
% 
%     folder_path = strcat(path, "\", foldername);
%  
%     [freq, data_raw] = import_data_4_fun(folder_path);
%     data(i,:) = data_raw(end,:);
% end
load NOCKIT5_2traces_data.mat
% data_mag = (abs(data).^2);
% data_dB = 10*log10(data_mag);


if smoothing
    for i=1:4
   data_dB(i,:) = smooth(data_dB(i,:), smoothing); 
    end
end

 % down sampling
res = 201;
%data_dB_red = nan(2,res);
freq_red = downsample(freq,res);
data_dB_red = (downsample(data_dB, res))';
% for i=1:2
%     
%  data_dB_red(:,i) = downsample(data_dB(:,i), res);
% end

%% plot 
figure(301); clf; plot(freq_red, data_dB_red(:,:), 'linewidth', 1.5); grid on; xlabel('freq (Hz)'); ylabel('dB'); title('measured transmission')
legend("4-->1","4-->2","4-->3","4-->4");
%  c_o = colororder;
 
%  colororder('default')
%% fitting
    v_ph = 1.8628e+06;
    v_ph_c = 5.4916e+06;
    Z0 = 50.1153;
    Zc = 1.9916e+03;

Y0 = 1/Z0;
Yc = 1/Zc;


% geometry: and network structure
N=30; % number of couplers. (= number of unit cells minus 1) 
M = 2; % number of lines
L0 = 100e-6; % length of each line segment
d = 26e-6; % length of each coupler segment

W=2.3e-6; % width of primary and secondary transmission lines
t=8e-9; % thickness of WSi (sputtered)
H=16e-9; % height of dielectric (say, Si - evaporated)
W_c=280e-9; % width of coupling line


[Y0, v_ph] = get_microstrip_properties(W,t,H);
[Yc, v_ph_c] = get_microstrip_properties(W_c,t,H);


% frequency etc.
%freq = 1e9*linspace(3,9,201); 
omega= 2*pi*freq;



input_idx = [1];   % can be more than one.
%%
% define graph: define an array of nodes with M rows and N+2 columns. the
% nodes are numbered such that nodes 1:M are the first column, M+1:2*M are
% the socond column  etc.
nodes = reshape(1:M*(N+2),M,N+2 );  
G = digraph();
% define main lines:
% connect each node (m,n) to (m, n+1), with weight  2.
% by convention weight 2 means a "main" edge and not a coupler.
s = nodes(:,1:N+1);
t = nodes(:,2:N+2);
w = 2*ones(size(s));
G = G.addedge(s,t,w);

% define coupling lines:
% connect each node (m,n) to (n+1,n) with weight 1.
% by convention weight 1 means a coupling edge.
s = reshape(nodes(1:M-1, 2:N+1), 1,(N)*(M-1));
t = reshape(nodes(2:M, 2:N+1), 1,(N)*(M-1));
w = ones(size(s));
G = G.addedge(s,t,w);


edge_num = G.numedges;
nodes_num = G.numnodes;

% figure(504)
% clf
% G.plot('xdata', x, 'ydata',y, 'linewidth', LWidths);

%


% define edges attributres; pahe velocity, length and characteristic
% admittance:
% rememeber weight = 1 means coupler edge, weight = 2 means regular edge:
%clearvars G.Edges.v_ph G.Edges.L G.Edges.Y
G.Edges.v_ph(G.Edges.Weight==2) = v_ph*ones(sum(G.Edges.Weight==2),1);
G.Edges.v_ph(G.Edges.Weight==1) = v_ph_c*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.Y(G.Edges.Weight==2) = Y0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Y(G.Edges.Weight==1) = Yc*ones(sum(G.Edges.Weight==1),1);

% define boundary conditions attribute for each edge according to the
% following convention:
% 0 - do nothing
% 1 - set t to 0
% 2 - set r to 0
% 3 - set t to 1
% 4 - set r t0 1
G.Edges.BC = zeros(edge_num,1);
G.Edges.BC(G.findedge(nodes(:,1),nodes(:,2))) = 1*ones(M,1);
G.Edges.BC(G.findedge(nodes(input_idx,1),nodes(input_idx,2))) = 3;
G.Edges.BC(G.findedge(nodes(:,N+1),nodes(:,N+2))) = 2*ones(M,1);


%%

graph_data = process_graph(G);






%% config
% options = optimoptions('fmincon','Display','iter', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');

options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','sqp-legacy', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','interior-point', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');
% options = optimoptions('fmincon','Display','iter','Algorithm','active-set', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');

% options = optimset('PlotFcns',@optimplotfval, 'Display','iter');
% %    x0 =  [1.99 ,   0.9 ,  0.8, -52];
% x0 = [1.9935    1.9193    0.8418  -21.9094];
% % x0 = [1.9, -52]
% lb = [.2,-60];
% ub = [3,-45];
% x0 = [2,2,1,-22]
x0 = [1,.8,1,1,1,-22];


 x0 =[1.0819    1.3596    0.4855    1.2106    1.0862  -22.0017]
 x0 = [0.9717    1.3580    0.4355    1.1422    1.0249]
 
    x0 = [1.0170    0.7107    0.2524    1.0619    0.9896]
%  x0 = [1,.8,.25]
 
%   x0 = [1.0163    0.7106    0.2522    1.0612    0.9889];
% x0 = [1,.8 , 0.25, 1,1];

x0 = [1,.8 , 0.25, -22];
x0= [  1.1051    0.7540    0.3074 ];

%   lb = [.6, .6, .1,0.8,0.7];
  lb = [.6, .6, .1];
% ub = [1.3, 1.4,2,1.2,1.3];
ub = [1.3, 1.4,2];

A = [];
b = [];
Aeq = [];
beq = [];

%% minimize
fun = @(x) get_cost_5(G,freq_red, data_dB_red, x);
fun(x0)
%%
[X,costval] = fmincon(fun,x0, A,b,Aeq,beq,lb,ub,[], options)
%  [X,costval] = fminsearchbnd(fun,x0,lb,ub ,options)

%% plot
% X = [1.6973    2.9982    2.9024  -52.7307];
% X = [2.1635    0.7438    0.7321  -51.2686];
% X = [2.1740    0.4795    0.4708,-51.2686]
% X = [2.1386    2.8675    2.9995];
% X = [1.0096    0.2403    0.2425];
%  X = [1.9935    0.9193    0.8418  -21.9094];
% offset = -51.2686;
figure(333); clf;
% 
%      v_ph_new = b_ph*X(1);
%     v_ph_c+_new = v_ph_c*X(2);
% %     Z0 = 50.1153;
%      Yc_new = Yc/X(3);

trans_dB = get_trans(G,freq,X, -22, 1:2);
%%
figure(333)
plot(freq, (trans_dB), 'linewidth',1.5 ,'linestyle', '--');
% ylim('auto')
% legend("4-->1","4-->2","4-->3","4-->4");
 grid on; xlabel('freq (Hz)'); ylabel('dB'); title('simulated transmission')
 hold on
  cc = colororder();
    plot(freq, data_dB(:,1), 'color',cc(1,:))
    plot(freq, data_dB(:,2), 'color',cc(2,:))
%  c_o = colororder;
 
%  colororder(c_o)
 


