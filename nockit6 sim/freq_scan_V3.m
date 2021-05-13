clearvars
tic 
% written by Guy 2020_08_31
% version 3 includes also attenuated reflections and loss
% edit 2020_11_15:
% parameters for nockit 6 design
% frequency scan 
Frequency=linspace(2,9,801)*1e9;
w = 2*pi*Frequency;

%% Configurations
coplanar_couplers =true; % toggle wether to calculate the coupler capacitance as a microstrip or as a coplanar
use_nockit_5_fit = false; 
% Lines gerometry:
    W=3e-6; % width of primary and secondary transmission lines
    t=8.5e-9; % thickness of WSi (sputtered)
    H=35e-9; % height of dielectric (say, Si - evaporated)
    W_c=200e-9; % width of coupling line
    gap_c = 5.4e-6; % coplanar gap for coupler. only relevant for the coplanar couplers case
 
% network geometry
    L = 100e-6; % length of each unit cell along main lines (m)
    d = 27e-6; % length of each coupling segment (m)
    N=32; % number of unit cells (31 couplers)
    M = 7; % number of lines
    idx_of_input_lines =4; % either an integer in [1,M]  or an array of integers in [1,M]

    if any(idx_of_input_lines > M)
        error("idx_of_imput_line greater than number of lines")
    end

% Electromagnetic properties:
    eps_r=11.7; % relative dielectric constant for Si
    eps_0=8.85e-12;

    % The kinetic inductance per unit length (L_kin) is calibrated according to measurement from 11.2.19 of a 10 nm / 2 micron strip
    L_kin=30.75615e-6*2e-6/W*10e-9/t;
    L_kin_c=30.75615e-6*2e-6/W_c*10e-9/t; % for coupling line
    
    % Formula for geometric inductance per unit length taken from https://www.allaboutcircuits.com/tools/microstrip-inductance-calculator/
    % (note that one must convert from inches). But L_geo<L_kin, so it might not be so important...
    L_geo=0.00508*39.3701*(log(2/(W+H))+0.5+0.2235*(W+H))*0.000001;
    

    C=W*eps_0*eps_r/H; % capacitance per unit length
    if coplanar_couplers
        % copied from CPWR_calculations.m:
        c = 3e8;
        u0 = 4*pi*1e-7;
        e_eff = (eps_r+1)/2*1.03^2;
       % Elliptic integrals
        k1 = W_c./(W_c+2*gap_c);
        k2 = sqrt(1-k1.^2);
        K1 = ellipke(k1.^2);
        K2 = ellipke(k2.^2);
        % Inductance and capacitance per unit length
        L_geo_c = u0/4*K2./K1;
        C_c = 4/(u0*c^2)*e_eff*K1./K2;
 
    else
        C_c=W_c*eps_0*eps_r/H; % capacitance per unit length for coupling line
        L_geo_c=0.00508*39.3701*(log(2/(W_c+H))+0.5+0.2235*(W_c+H))*0.000001;
    end

    L_tot=L_geo+L_kin; % total inductance per unit length;
    L_tot_c=L_geo_c+L_kin_c; % total inductance per unit length for coupling line;

% loss mechanisms:
    R = 0e4; % series resistance per unit length  
    G = .5e1; % shunt conductance per unit length
    R_c =0;  % series resistance per unit length for couplers 
    G_c =0;  % shunt conductance per unit length for couplers

% Phase velocity and characteristic impedance
    v_ph=1/sqrt(L_tot*C);
    v_ph_c=1/sqrt(L_tot_c*C_c);
    Z_0=sqrt(L_tot/C);
    Z_c=sqrt(L_tot_c/C_c);

    
    
 if use_nockit_5_fit
% parameters correction from fit. use these to get somthing close to the
% measurement for 2 traces NOCKIT5, but note that we still have to explain the factor of 2 in
% the phase velocity. the other two factors are close to 1, so they are OK.
x = [1.9935    0.9193    0.8418];     
%x = [2,1,1]
v_ph = v_ph*x(1);
    v_ph_c = v_ph_c*x(2);
    Z_c =  Z_c*x(3);
    
 end

% boundaries reflections and attenuators:
%  see pdf for details
    
    % attenuation in dB at the interfaces:  first element is for input
    % ports, second element is for output ports. set to 0 for full
    % reflections, set to inf for no reflections
        atten_dB =[inf,inf] ; 
        
    atten = 10.^(-atten_dB/20); % voltage attenuation factor
    coax_l =30e-2; % estimated total length of coax cable. set to 0 for reflections from the chip itself
    v_ph_coax = (2/3)*3e8; % phase velocity in coax lines

% complex wavenumbers and and impedances

k=1i*sqrt((R-1i*w*Z_0/v_ph).*(G-1i*w/(Z_0*v_ph))); % effective complex wave number for primary lines (rad/m)
k = k.*sign(real(k)); % becuase the sqrt() makes a choice and we want the one where real part of k to be positive. see pdf.
k_c=1i*sqrt((R_c-1i*w*Z_c/v_ph_c).*(G_c-1i*w/(Z_c*v_ph_c))); % effective complex wave number for coupling lines (rad/m)
k_c= k_c.*sign(real(k_c)); % becuase the sqrt() makes a choice and we want the one where real part of k to be positive. see pdf.
Zeff = sqrt((R-1i*w*Z_0/v_ph)./(G-1i*w/(Z_0*v_ph)));
Zeff_c = sqrt((R_c-1i*w*Z_c/v_ph_c)./(G_c-1i*w/(Z_c*v_ph_c)));

k_coax = 2*pi*Frequency/v_ph_coax;  % wave number for coax lines (rad/m)

Y_0 = 1./Zeff;
Y_c = 1./Zeff_c;

%% loop on frequencies

% pre-allocate
    trans  = nan(length(Frequency), M);
    ref  = nan(length(Frequency), M);

%loop
for idx = 1:length(Frequency)
    
    ref_factor = (atten.^2).*exp(2i*k_coax(idx)*coax_l); % a complex factor that will be scaling the reflections at the boundaries

    
    [t, r] = SolveArrayFunV3(M,N, L,d,k(idx),k_c(idx),Y_0(idx), Y_c(idx),idx_of_input_lines,ref_factor);
    trans(idx,:) = t(:,end);
    ref(idx,:) = r(:,1);
    
end

trans_mag = abs(trans);
ref_mag = abs(ref);

trans_phase = unwrap(angle(trans));

trans_dB = 20*log10(trans_mag);
ref_dB = 20*log10(ref_mag);

toc
%% plot

figure(1) % transmittance graphs
    plot(Frequency, trans_dB(:,1:4)', 'linewidth', 1.5)
    grid on
    leg = legend(num2str((1:4)'),"location", "best", "fontsize", 16);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("dB", "fontsize", 16)
    title_str = sprintf('lines transmission w. input from line %d \n NOCKIT5 fit = %d \n coplanar couplers = %d' , idx_of_input_lines, use_nockit_5_fit,coplanar_couplers);
    title(title_str, "fontsize", 16)

figure(2)  % reflectance graphs  
    plot(Frequency, ref_dB(:,1:4)', 'linewidth', 1.5)
    grid on
    leg = legend(num2str((1:4)'),"location", "best", "fontsize", 16);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("dB", "fontsize", 16)
    title("reflectance", "fontsize", 16)

figure(3) % colormaps
   subplot(2,1,1) % transmittance
       imagesc( trans_dB', "XData", Frequency, "YData", (1:M))
       view(2)
       shading flat
       colorbar
       xlabel("frequency (Hz)", "fontsize", 13)
       ylabel("line", "fontsize", 16)
       title("transmittace (dB)", "fontsize", 16)  
    subplot(2,1,2) % reflectance
       imagesc( ref_dB', "XData", Frequency, "YData", (1:M))
       view(2)
       shading flat
       colorbar
       xlabel("frequency (Hz)", "fontsize", 13)
       ylabel("line", "fontsize", 16)
       title("reflectance (dB)", "fontsize", 16)
    
       colormap jet
%%

figure(4)
clf
    plot(Frequency, trans_phase(:,4)', 'linewidth', 1.5)
    grid on
    leg = legend(num2str((4)'),"location", "best", "fontsize", 16);
    title(leg, "line")
    xlabel("frequency (Hz)", "fontsize", 16)
    ylabel("phase", "fontsize", 16)
    title_str = sprintf('phase at output w. input from line %d \n NOCKIT5 fit = %d \n coplanar couplers = %d' , idx_of_input_lines, use_nockit_5_fit,coplanar_couplers);
    title(title_str, "fontsize", 16)

    
 