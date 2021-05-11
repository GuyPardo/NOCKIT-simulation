clearvars




%% import data
smoothing = 0;

traceid = 1:4;
path =  "C:\Users\guypa\Google Drive\LIMUDIM\Lab_research\nockit6 data\Nockit6\Transmission@-40";
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
figure(1); clf; plot(freq_red, data_dB_red(:,:), 'linewidth', 1.5); grid on; xlabel('freq (Hz)'); ylabel('dB'); title('measured transmission')
legend("4-->1","4-->2","4-->3","4-->4");
 c_o = colororder;
 
 colororder('default')
%% fitting

%% config
options = optimset('PlotFcns',@optimplotfval, 'Display','iter');

  x0 =  [1.99 ,   0.9 ,  0.8, -52];



  lb = [.2, .2, .2,-70];
ub = [3, 3,3,-40];

A = [];
b = [];
Aeq = [];
beq = [];

%% minimize
fun = @(x) get_cost_4(freq_red, data_dB_red, x);
[X,costval] = fmincon(fun,x0, A,b,Aeq,beq,lb,ub,[] ,options)

%% plot
X = [1.6973    2.9982    2.9024  -52.7307];

figure(2); clf;

    v_ph = 1.8628e+06*X(1);
    v_ph_c = 5.4916e+06*X(2);
    Z0 = 50.1153;
    Zc = 1.9916e+03*X(3);

trans = freq_scan_fun_4(freq,v_ph, v_ph_c, 1/Z0, 1/Zc );
%%

plot(freq, 20*log10(abs(trans(4,:))) + X(4), 'linewidth',1.5 );
% ylim('auto')
legend("4-->1","4-->2","4-->3","4-->4");
 grid on; xlabel('freq (Hz)'); ylabel('dB'); title('simulated transmission')
%  c_o = colororder;
 
 colororder(c_o)
 


