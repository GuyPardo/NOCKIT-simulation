clearvars
addpath('')
name = "\Transmission@-40";

lab = true;


if lab

   folder = 'Z:\Measurements\2020.dir\01_12_20 cooldown\Nockit6';
else
   folder = 'C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\nockit6 data\Nockit6'; 
end
%% import data
smoothing = 0;

traceid = 1:4;
path =  strcat(folder,name);
clearvars data
for i=1:length(traceid)

        foldername = sprintf("4to%s", num2str(traceid(i)));

    folder_path = strcat(path, "\", foldername);
 
    [freq, data_raw] = import_data_4_fun(folder_path);
    data(i,:) = data_raw(end,:);
end

data_mag = (abs(data).^2);
data_dB = 10*log10(data_mag);


if smoothing
    for i=1:4
   data_dB(i,:) = smooth(data_dB(i,:), smoothing); 
    end
end

 % down sampling
res = 201;
data_dB_red = nan(4,res);
freq_red = downsample(freq,res);
for i=1:4
    
 data_dB_red(i,:) = downsample(data_dB(i,:), res);
end

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
N=31; % number of couplers. (= number of unit cells minus 1) 
M = 7; % number of lines
L0 = 100e-6; % length of each line segment
d = 27e-6; % length of each coupler segment



% frequency etc.
%freq = 1e9*linspace(3,9,201); 
omega= 2*pi*freq;



input_idx = [4];   % can be more than one.
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
options = optimoptions('fmincon','Display','iter','Algorithm','sqp', 'PlotFcn', 'optimplotfval','PlotFcn', 'optimplotfval');

   x0 =  [1.99 ,   0.9 ,  0.8, -52];
% % x0 = [1.9, -52]
% lb = [.2,-60];
% ub = [3,-45];

  lb = [.2, .2, .2,-60];
ub = [3, 3,3,-45];

A = [];
b = [];
Aeq = [];
beq = [];

%% minimize
fun = @(x) get_cost_5(G,freq_red, data_dB_red, x);
[X,costval] = fmincon(fun,x0, A,b,Aeq,beq,lb,ub,[], options)

%% plot
% X = [1.6973    2.9982    2.9024  -52.7307];
X = [2.1635    0.7438    0.7321  -51.2686];
% X = [2.1740    0.4795    0.4708,-51.2686]
% X = [2.1386    2.8675    2.9995];
X = [1.0096    0.2403    0.2425];

offset = -51.2686;
figure(333); clf;
% 
%     v_ph = 1.8628e+06*X(1);
%     v_ph_c = 5.4916e+06*X(2);
%     Z0 = 50.1153;
%      Zc = 1.9916e+03*X(3);

trans_dB = get_trans(G,freq,X, offset, 1:4);
%%
figure(333)
plot(freq, (trans_dB(1:4,:)), 'linewidth',1.5 );
% ylim('auto')
legend("4-->1","4-->2","4-->3","4-->4");
 grid on; xlabel('freq (Hz)'); ylabel('dB'); title('simulated transmission')
%  c_o = colororder;
 
%  colororder(c_o)
 


