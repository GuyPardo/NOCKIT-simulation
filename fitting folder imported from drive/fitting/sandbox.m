% import data

path = 'Z:\Measurements\2020.dir\13_07_20 cooldown\Mosquito\Network analyzer';

[freq, dir] = ReadCTI(path, 'in4_out2');
[freq, cop] = ReadCTI(path, 'in4_out1');

dir_dB = 20*log10(abs(dir));
cop_dB = 20*log10(abs(cop));

%%
%default  config
N = 31;
L = 100e-6;
d = 20e-6;

% default config
v_ph  = 1.6850e+06;
v_ph_c = 1.7626e+06;
Z_0 = 49.8391;
Z_c = 547.9262;

freqArr = freq;

% % scanning changes in the default config
% factor = linspace(0.5,2, 10);
% NN = length(factor);
% 
% cost = nan(1,NN);
% for i = 1:NN
% 
%     
%                 cost(i) = get_cost(factor(i));
% 
% end
%   %%
%% down sampling


idx = 1:length(freq);                                 % Index
idxq = linspace(min(idx), max(idx), 201);    % Interpolation Vector
freq_red = interp1(idx, freq, idxq, 'linear');       % Downsampled Vector
dir_dB_red = interp1(idx, dir_dB, idxq, 'linear');
cop_dB_red = interp1(idx, cop_dB, idxq, 'linear');


options = optimset('PlotFcns',@optimplotfval, 'Display','iter');
x0 = [1.2,1.2,0.8,-15];
% x0 =  [1.6773    0.7252    0.5682  -21.8800,0,0];
  x0 = [1.6963    0.7016    0.5538  0.01]
  x0 =  [1.9868    0.9295    0.8616, -20];

  lb = [.4, .4, .4,-50];
ub = [2, 2, 2,0];
%   x0 = [1.2  ,  1.2   ,0.8  0.01,0,0];
%% matlab optimization


%   [X,costval] = fminsearchbnd(@get_cost,x0,lb,ub, options)
% 
A = [];
b = [];
Aeq = [];
beq = [];



[X,costval] = fmincon(@get_cost_3,x0, A,b,Aeq,beq,lb,ub,[] ,options);



%% fitting
lb = [.4, .4, .4, -30];
ub = [2,2,2, 0];
% X = lsqcurvefit(@coupled, x0,freq_red, cop_dB_red, lb,ub, optimoptions('lsqcurvefit','Algorithm','trust-region-reflective','PlotFcns',@optimplotfval, 'Display','iter'))

%%


figure
plot(factor,cost)
factor(3)


%%
% X = [1.2047    1.2003    0.8004  -21.7684]
% X = [1.6963    0.7016    0.5538  -21.8633]
% X = [1.6417    1.9999    0.9151]
freq2 = linspace(0,10,2001)*1e9;
[direct_trans,coupled_trans] = freq_scan_fun(freq2,N,L,d,v_ph*X(1),v_ph_c*X(2),Z_0,Z_c*X(3));
 close all
figure(22)
plot(freq, (dir_dB))
hold on
plot(freq, (cop_dB))
legend(["direct", "coupled"], 'fontsize', 16)
ylabel("dB", 'fontsize', 16)
xlabel ('frequency (Hz)', 'fontsize', 16 )
title ('Measured Transmission', 'fontsize', 16 )

grid on
figure(23)
plot(freq2, 20*log10(abs(direct_trans)))
hold on
plot(freq2, 20*log10(abs(coupled_trans)))
grid on
legend(["direct", "coupled"], 'fontsize', 16)
ylabel("dB", 'fontsize', 16)
xlabel ('frequency (Hz)', 'fontsize', 16 )
title ('Simulated Transmission', 'fontsize', 16 )

%%
close all
 [pks, locs ] = findpeaks(- 20*log10(abs(direct_trans)), 'MinPeakProminence',10) 
 findpeaks(- 20*log10(abs(coupled_trans)), 'MinPeakProminence',10) 
