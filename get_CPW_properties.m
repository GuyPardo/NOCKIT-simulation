function [Y0, v_ph] = get_CPW_properties(thickness, width,gap)
% returns the transmission line properties (admittance and phase-velocity)
% for a WSi coplanar wave guide transmission line

% Lines gerometry:
    W=width; 
    t=thickness; 



% Electromagnetic properties:
    eps_r=11.7; % relative dielectric constant for Si
    eps_0=8.85e-12;

% The kinetic inductance per unit length (L_kin) is calibrated according to measurement from 11.2.19 of a 10 nm / 2 micron strip
    L_kin=30.75615e-6*2e-6/W*10e-9/t;


        % copied from CPWR_calculations.m:
        c = 3e8;
        u0 = 4*pi*1e-7;
        e_eff = (eps_r+1)/2*1.03^2;
       % Elliptic integrals
        k1 = W./(W+2*gap);
        k2 = sqrt(1-k1.^2);
        K1 = ellipke(k1.^2);
        K2 = ellipke(k2.^2);
        % Inductance and capacitance per unit length
        L_geo = u0/4*K2./K1;
        C = 4/(u0*c^2)*e_eff*K1./K2;


    L_tot=L_geo+L_kin; % total inductance per unit length;


    
    

% Phase velocity and characteristic impedance
    v_ph=1/sqrt(L_tot*C);

    Z_0=sqrt(L_tot/C);

    


 
Y0 = 1/Z_0;

 
 
end
