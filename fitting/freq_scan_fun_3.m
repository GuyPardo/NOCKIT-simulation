function [trans] = freq_scan_fun_3(FrequencyArr,M,N,L,d,v_ph,v_ph_c,Y_0,Y_c )
%UNTITLED4 Summary of this function goes here

kArr=2*pi*FrequencyArr/v_ph; % wave number for primary lines (rad/m)
k_cArr=2*pi*FrequencyArr/v_ph_c; % wave number for coupling segmets (rad/m)

% pre allocation:
% t_p_finalArr = nan(size(FrequencyArr));
% t_s_finalArr = nan(size(FrequencyArr));
% r_p_1Arr = nan(size(FrequencyArr));
% r_s_1Arr = nan(size(FrequencyArr));

trans = nan(M,length(FrequencyArr));

% loop
for idx = 1:length(FrequencyArr)
    k = kArr(idx);
    k_c = k_cArr(idx);
    
    [t,r] = SolveArrayFunV2(M,N, L,d,k,k_c,Y_0, Y_c,1);
    
    trans(:,idx) = t(:,end);

    
end



end

