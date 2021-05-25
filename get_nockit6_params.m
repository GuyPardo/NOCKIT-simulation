function nockit_params = get_nockit6_params(X)

% nockit6  params

nockit_params.N=31; % number of couplers. (= number of unit cells minus 1) 
nockit_params.M = 7; % number of lines
nockit_params.L0 = 100e-6; % length of each line segment
nockit_params.d = 27e-6; % length of each coupler segment

nockit_params.W=3e-6; % width of primary and secondary transmission lines
nockit_params.t=8.5e-9; % thickness of WSi (sputtered)
nockit_params.H=29e-9; % height of dielectric (say, Si - evaporated)
nockit_params.W_c=200e-9; % width of coupling line
nockit_params.gap_c = 8e-6;
nockit_params.gnd_cond=0;

if nargin>0
    nockit_params.t = nockit_params.t*X(1);
    nockit_params.W = nockit_params.W*X(2);
    nockit_params.W_c = nockit_params.W_c*X(3);
    nockit_params.H = nockit_params.H*X(4);
    
end
end