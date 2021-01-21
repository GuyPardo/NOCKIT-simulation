function [down_sampled_vec] = downsample(vec, new_res)
%UNTITLED9 Summary of this function goes here
%   Detailed explanation goes here
idx = 1:length(vec);                                 % Index
idxq = linspace(min(idx), max(idx), new_res);    % Interpolation Vector
down_sampled_vec = interp1(idx, vec, idxq, 'linear');       % Downsampled Vector

end

