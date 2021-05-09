% fiting nockit 2

% load('C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\NOCKIT_new_folder\data\NOCKIT5_2traces_data.mat')
load NOCKIT5_2traces_data
N_downsample = 201;
%N_downsample = length(freq);
data_red = downsample(data,N_downsample);
data_dB_red = 20*log10(abs(data_red));
freq_red = downsample(freq,N_downsample);
figure(1); plot(freq_red, data_dB_red); title("downlampled data"); grid on;

dir = data_dB_red(:,1);
cop = data_dB_red(:,2);
%  lcs = find(abs(diff(cop))>1);

%  lcs(1)

 findpeaks(-data_dB_red(:,1), 'MinPeakProminence',10)
%  findpeaks(-data_dB_red(:,2), 'MinPeakProminence',10)
%%
sig_pwr = -80; %dbm
% geometry: and network structure
N=30; % number of couplers. (= number of unit cells minus 1) 
M = 2; % number of lines
L0 = 100e-6; % length of each line segment
d = 20e-6; % length of each coupler segment
t = 9e-9;
W = 2.3e-6;
W_c = 300e-9;
H = 16e-9;



addpath('C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\repos2\NOCKIT-simulation');
addpath('Z:\Users\Guy\coupling transission lines\repos\NOCKIT-simulation');
[Y0, v_ph]  = get_microstrip_properties(W,t,H);
[Yc, v_ph_c]  = get_microstrip_properties(W_c,t,H);

%%

x0 = [1.7,1.7,1.1,-20];
% x0 = [1.7,1.7,1.1]
% x0 = [1.7,1.1];
% x0 = [1.7,1.1, -22];



x0 = [1,1,2,1];

%  x0 = [0.9490    1.1236    2.6439  -22.2561]



% x0=   [0.7938    0.8354, 1    1.9709]
% x0 = [1,1,2]
x0 = [0.7938    0.8354    1.9709, -22  ]

x0 = [1.9,1.1,1/0.8, -22]
 x0 =  [1.9935    0.9193    1/0.8418  -21.9094]
% 
% x0 = [0.7948    0.8416    1.9652    1.1838]
% x0 = [0.7948    0.8416    1.9652, 1*20/27]

%  x0 = [0.7552    0.9859    2.0450    1.2034]
%  x0 = [  0.7898    1.1247    2.0031    1.4134]
%x0 = [1,1,1]
cost = get_cost_non_lin(x0)

%%
%% matlab optimization


options = optimset('PlotFcns',@optimplotfval, 'Display','iter');

%   lb = [.5, .5, .1 .7];
  lb = [0.5 0.5, 0.5 ,-35 ];
%    lb = [.1, .1, -40];
% ub = [2.5, 2 ,3,1.3];
ub = [3, 3,3, 0];
% ub = [3, 3, 0];
%     [X,costval] = fminsearchbnd(@get_cost_non_lin,x0,lb,ub, options)

A = [];
b = [];
Aeq = [];
beq = [];



 [X,costval] = fmincon(@get_cost_non_lin,x0, A,b,Aeq,beq,lb,ub,[] ,options);
% % 

%%
% X = [ 0.1216    0.9369    1.0754  -22.0307]
%X = [1.7,1.1]
% X=[1.7856    1.1320  -22.0302]
% 
%  X = [1 ,1,   2  -22.2021]
%  X=[0.794    0.8454    1.9609, 1];
%     X = [.7938    0.8354    1.9709]
% X =     [0.7946    0.8689    1.8981 ]
% X =  [1.2418    0.7428    1.0084]
%  X = [  0.7864    0.8607    1.9643]
 
%  X = [0.8    0.8007    1.943]
% X = [0.7252    0.9959    2.0450    1.3034]
% X = [ 0.7898    1.1247    2.0031    1.4134]
% X = [0.7600    0.9871    2.0499    1.2147]
%  X = [0.7668    1.0    1.91    1.1406];
X = [ 1.9935    0.9193    1/0.8418  -21.9094]
t_new = t;%*X(1);

W_c_new = W_c;%*X(2);
W_new = W;%*X(3);

d_new=d;%*X(4);
[Y0_new, v_ph_new]  = get_microstrip_properties(W_new,t_new,H);
[Yc_new, v_ph_c_new]  = get_microstrip_properties(W_c_new,t_new,H);
%     v_ph_new = v_ph_new*X(3);
%     v_ph_c_new = v_ph_c_new*X(3);
%     Yc_new = Yc_new*X(3);
%     Y0_new = Y0_new*X(3);
    
    
       v_ph_new = v_ph_new*X(1);
    v_ph_c_new = v_ph_c_new*X(2);
    Yc_new = Yc_new*X(3);
    Y0_new = Y0_new;%*X(3);
    
      
 Ic = 60e-6*(t_new/6e-9)*(W/2.3e-6);%*X(4); % % from mikita measurement
Icc = 60e-6*(t_new/6e-9)*(W_c_new/2.3e-6);%*X(4); % from mikita measurement
    trans = freq_scan_non_lin_fun(freq_red,v_ph_new, v_ph_c_new, Y0_new, Yc_new,Ic,Icc,sig_pwr,d_new);
    trans_dB = 20*log10(abs(trans));
%
    figure(2); clf;plot(freq_red, trans_dB - 21 , '--') ; grid on
    cc = colororder();
    hold on
    plot(freq, data_dB(:,1), 'color',cc(1,:))
    plot(freq, data_dB(:,2), 'color',cc(2,:))
%     findpeaks(-trans_dB(2,:), 'MinPeakProminence',3);
% [peaks2_sim, locs2_sim] = findpeaks(-data_dB_red(2,:), 'MinPeakProminence',10);
dir = trans_dB(1,:);
cop = trans_dB(2,:);
%%
BG1 = freq_red(find(abs(diff(dir))>2,1));
BG2 = freq_red(find(abs(diff(cop))>2,1));
    