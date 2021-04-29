% fiting nockit 2

%load('C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\NOCKIT_new_folder\data\NOCKIT5_2traces_data.mat')
load NOCKIT5_2traces_data
N_downsample = 201;
N_downsample = length(freq);
data_red = downsample(data,N_downsample);
data_dB_red = 20*log10(abs(data_red));
freq_red = downsample(freq,N_downsample);
figure(1); plot(freq_red, data_dB_red); title("downlampled data"); grid on;


%%
% geometry: and network structure
N=30; % number of couplers. (= number of unit cells minus 1) 
M = 2; % number of lines
L0 = 100e-6; % length of each line segment
d = 27e-6; % length of each coupler segment
t = 9e-9;
W = 2.3e-6;
W_c = 300e-9;
H = 16e-9;


addpath('C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\repos2\NOCKIT-simulation');
[Y0, v_ph]  = get_microstrip_properties(W,t,H);
[Yc, v_ph_c]  = get_microstrip_properties(W_c,t,H);

%%

x0 = [1.7,1.7,1.1,-20];
% x0 = [1.7,1.7,1.1]
% x0 = [1.7,1.1];
% x0 = [1.7,1.1, -22];



x0 = [1,1,2,-20];

x0 = [0.9490    1.1236    2.6439  -22.2561]
cost = get_cost_new(x0)

%%
%% matlab optimization


options = optimset('PlotFcns',@optimplotfval, 'Display','iter');

  lb = [.5, .5, .1 ,-40 ];
%   lb = [.1, .1, .1  ];
%    lb = [.1, .1, -40];
ub = [2.5, 2,3, 0];
% ub = [3, 3,3];
% ub = [3, 3, 0];
%    [X,costval] = fminsearchbnd(@get_cost_new2,x0,lb,ub, options)
% 
A = [];
b = [];
Aeq = [];
beq = [];



[X,costval] = fmincon(@get_cost_new2,x0, A,b,Aeq,beq,lb,ub,[] ,options);


%%
% X = [ 0.1216    0.9369    1.0754  -22.0307]
%X = [1.7,1.1]
% X=[1.7856    1.1320  -22.0302]
% 
%  X = [1 ,1,   2  -22.2021]

t_new = t*X(1);
W_c_new = W_c*X(2);
[Y0_new, v_ph_new]  = get_microstrip_properties(W,t_new,H);
[Yc_new, v_ph_c_new]  = get_microstrip_properties(W_c_new,t_new,H);
    v_ph_new = v_ph_new*X(3);
    v_ph_c_new = v_ph_c_new*X(3);
    Yc_new = Yc_new*X(3);
    Y0_new = Y0_new*X(3);
    
      
    
    trans = freq_scan_fun_4(freq_red,v_ph_new, v_ph_c_new, Y0_new, Yc_new);
    trans_dB = 20*log10(abs(trans));
%
    figure(2); clf;plot(freq_red, trans_dB + X(4), '--') ; grid on
    cc = colororder();
    hold on
    plot(freq, data_dB(:,1), 'color',cc(1,:))
    plot(freq, data_dB(:,2), 'color',cc(2,:))
    
    