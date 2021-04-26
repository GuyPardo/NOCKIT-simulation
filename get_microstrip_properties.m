function [Y0, v_ph] = get_microstrip_properties(width,metal_thickness,dielectric_thickness)
% returns the transmission line properties (admittance and phase-velocity)
% for a WSi microstrip transmission line

% Lines gerometry:
    W=width; % width of primary and secondary transmission lines
    t=metal_thickness; % thickness of WSi (sputtered)
    H=dielectric_thickness; % height of dielectric (say, Si - evaporated)
    
    
 

% Electromagnetic properties:
    eps_r=11.7; % relative dielectric constant for Si
    eps_0=8.85e-12;

    % The kinetic inductance per unit length (L_kin) is calibrated according to measurement from 11.2.19 of a 10 nm / 2 micron strip
    L_kin=30.75615e-6*2e-6/W*10e-9/t;
  
    
    % Formula for geometric inductance per unit length taken from https://www.allaboutcircuits.com/tools/microstrip-inductance-calculator/
    % (note that one must convert from inches). But L_geo<L_kin, so it might not be so important...
    L_geo=0.00508*39.3701*(log(2/(W+H))+0.5+0.2235*(W+H))*0.000001;
%     L_geo = 0;

    C=W*eps_0*eps_r/H; % capacitance per unit length


    L_tot=L_geo+L_kin; % total inductance per unit length;
  

    
    

% Phase velocity and characteristic impedance
    v_ph=1/sqrt(L_tot*C);

    Z_0=sqrt(L_tot/C);

 
Y0 = 1/Z_0;

 
end
