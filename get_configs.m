function [Y0, v_ph] = get_CPW_properties(W,t,H,W_c,gap_c,coplanar_couplers, use_notckit_5_fit) 
% returns the transmission line properties from the geometric 

% Lines gerometry:
    W=3e-6; % width of primary and secondary transmission lines
    t=10e-9; % thickness of WSi (sputtered)
    H=16e-9; % height of dielectric (say, Si - evaporated)
    W_c=200e-9; % width of coupling line
    gap_c = 1.4e-6; % coplanar gap for coupler. only relevant for the coplanar couplers case
 

% Electromagnetic properties:
    eps_r=11.7; % relative dielectric constant for Si
    eps_0=8.85e-12;

    % The kinetic inductance per unit length (L_kin) is calibrated according to measurement from 11.2.19 of a 10 nm / 2 micron strip
    L_kin=30.75615e-6*2e-6/W*10e-9/t;
  


        % copied from CPWR_calculations.m:
        c = 3e8;
        u0 = 4*pi*1e-7;
        e_eff = (eps_0+1)/2*1.03^2;
       % Elliptic integrals
        k1 = W_c./(W_c+2*gap_c);
        k2 = sqrt(1-k1.^2);
        K1 = ellipke(k1.^2);
        K2 = ellipke(k2.^2);
        % Inductance and capacitance per unit length
        L_geo = u0/4*K2./K1;
        C = 4/(u0*c^2)*e_eff*K1./K2;


    L_tot=L_geo+L_kin; % total inductance per unit length;
%     L_tot_c=L_geo_c+L_kin_c; % total inductance per unit length for coupling line;

    
    

% Phase velocity and characteristic impedance
    v_ph=1/sqrt(L_tot*C);

    Z_0=sqrt(L_tot/C);

    


 
Y0 = 1/Z_0;

 
 
end
