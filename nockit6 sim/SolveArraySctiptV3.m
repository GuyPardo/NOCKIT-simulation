%  written by Guy 2020_08_24
% This script solves the MxN network for a single frequency (6 GHz by default),
% reconstructs the voltage and current, and plots the power propagation
% along the lines.
% version 3 can also support loss and attenuated reflections

% for loops on parameters (like frequency) it is more convenient to use
% SolveArrayFunV3.m
clearvars
%% Configurations
coplanar_couplers =false; % toggle wether to calculate the coupler capacitance as a microstrip or as a coplanar
phase = transpose(linspace(0,2*pi,201));
phase_factor = [ones(size(phase)),exp(1i*phase)];

% Geometry:
% Lines gerometry:
W=3e-6; % width of primary and secondary transmission lines
t=10e-9; % thickness of WSi (sputtered)
H=16e-9; % height of dielectric (say, Si - evaporated)
W_c=200e-9; % width of coupling line
gap_c = 1.4e-6; 
% network geometry
L = 100e-6; % length of each unit cell along main lines (m)
d = 27e-6; % length of each coupling segment (m)
N=32; % number of unit cells
M =7; % number of lines
idx_of_input_lines = [4];
if any(idx_of_input_lines > M)
    error("idx_of_imput_line greater than number of lines")
end

% Electromagnetic properties:
eps_r=11.7; % relative dielectric constant for Si
eps_0=8.85e-12;

L_kin=30.75615e-6*2e-6/W*10e-9/t;
L_kin_c=30.75615e-6*2e-6/W_c*10e-9/t; % for coupling line
% The kinetic inductance per unit length (L_kin) is calibrated according to measurement from 11.2.19 of a 10 nm / 2 micron strip
L_geo=0.00508*39.3701*(log(2/(W+H))+0.5+0.2235*(W+H))*0.000001;
% Formula for geometric inductance per unit length taken from https://www.allaboutcircuits.com/tools/microstrip-inductance-calculator/
% (note that one must convert from inches). But L_geo<L_kin, so it might not be so important...
C=W*eps_0*eps_r/H; % capacitance per unit length



if coplanar_couplers
        % copied from CPWR_calculations.m:
        c = 3e8;
        u0 = 4*pi*1e-7;
        e_eff = (eps+1)/2*1.03^2;
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

    %loss mechanisms:
    R =0e4; % series resistance per unit length (Ohm/m)
    G = 0e1; % shunt conductance per unit length (1/Ohm*m)
    R_c =0; 
    G_c =0;


% Phase velocity, impedance and wave number:
v_ph=1/sqrt(L_tot*C);
v_ph_c=1/sqrt(L_tot_c*C_c);
Z_0=sqrt(L_tot/C);
Z_c=sqrt(L_tot_c/C_c);

v_ph=2*v_ph;

Y_0 = 1/Z_0;
Y_c = 1/Z_c;

Frequency=6e9; % THIS can be a parameter to run over later
k=2*pi*Frequency/v_ph; % wave number for primary lines (rad/m)
k_c=2*pi*Frequency/v_ph_c; % wave number for coupling segmets (rad/m)

wavelength = 2*pi/k;
wavelength_c = 2*pi/k_c;
% the wavelength variables are not used in the simulation, and are
% calculated here for convenience

% boundaries reflections and attenuators:
%  see pdf for detrails
atten_dB =[inf,inf]; % set to 0 for full reflections, and to inf for no reflections. component(1) is for inputs, componenent (2) is for output.

atten = 10.^(-atten_dB/20); % voltage attenuation factor
coax_l = 0e-2; % estimated total length of coax cable
v_ph_coax = (2/3)*3e8; % phase velocity in coax lines

k_coax = 2*pi*Frequency/v_ph_coax;
factor(1) = -(atten(1)^2)*exp(2i*k_coax*coax_l);
factor(2) = -(atten(2)^2)*exp(2i*k_coax*coax_l);
w = 2*pi*Frequency;
% effective complex propagation constant
k=1i*sqrt((R-1i*w*Z_0/v_ph).*(G-1i*w/(Z_0*v_ph))); % effective complex wave number for primary lines (rad/m)
k = k*sign(real(k));
k_c=1i*sqrt((R_c-1i*w*Z_c/v_ph_c).*(G_c-1i*w/(Z_c*v_ph_c)));
k_c = k_c*sign(real(k_c));

Zeff = sqrt((R-1i*w*Z_0/v_ph)./(G-1i*w/(Z_0*v_ph)));
Zeff_c = sqrt((R_c-1i*w*Z_c/v_ph_c)./(G_c-1i*w/(Z_c*v_ph_c)));


tic
%% encoding equations
%  see the pdf for details on the encoding procedure
eqn_num = 4*M*N - 2*M - 2*N +2;


%pre-allocate
big_ten = zeros(eqn_num + 2*M-2, 2*M - 1, 2*N); % need to check the size. 
eqn_count = 0; % count how many eqns we have encoded

% bulk: 
% encoding the 4 equations for a single bulk unit cell


unit_ten = nan(4,4,4);

unit_ten(1,:,:) = [0 0 0 0;
                   0 0 0 0;
                   0 0 -1 -1;
                   0 0 1 1];
unit_ten(2,:,:) = [0 0 0 0;
                   0 0 0 0;
                   exp(1i*k*L) exp(-1i*k*L) -1 -1; 
                   0 0 0 0 ];
unit_ten(3,:,:) = [0 0 0 0;
                   0 0 exp(1i*k_c*d) exp(-1i*k_c*d);
                   0 0 -1 -1;
                   0 0 0 0];
unit_ten(4,:,:) = [0 0 0 0;
                   0 0 -Y_c*exp(1i*k_c*d) Y_c*exp(-1i*k_c*d);
                   -Y_0*exp(1i*k*L) Y_0*exp(-1i*k*L) Y_0 -Y_0;
                   0 0 Y_c -Y_c];
               
               
% encoding the entire bulk equations

for j = 2:M-1
    for n = 2:N
        
        big_ten((1:4)+eqn_count  ,  (1:4) + (j-2)*2  ,  (1:4) + (n-2)*2) = unit_ten;
        eqn_count = eqn_count+4;             
    end
end


% j=1 and j=M  boundary: 
% encoding the 3 equations for a single cell in the j=1 boundary:
unit_ten_j_1 = nan(3,2,4);


unit_ten_j_1(1,:,:) = [0 0 -1 -1;
                   0 0 1 1];
unit_ten_j_1(2,:,:) = [exp(1i*k*L) exp(-1i*k*L) -1 -1; 
                   0 0 0 0 ];

unit_ten_j_1(3,:,:) = [-Y_0*exp(1i*k*L) Y_0*exp(-1i*k*L) Y_0 -Y_0;
                   0 0 Y_c -Y_c];
               
               
%  encoding the 3 equations for a single cell in the j=M boundary:
unit_ten_j_M = nan(3,3,4);


unit_ten_j_M(1,:,:) = [0 0 0 0;
                   0 0 0 0;
                   exp(1i*k*L) exp(-1i*k*L) -1 -1]; 

unit_ten_j_M(2,:,:) = [0 0 0 0;
                   0 0 exp(1i*k_c*d) exp(-1i*k_c*d);
                   0 0 -1 -1];

unit_ten_j_M(3,:,:) = [0 0 0 0;
                   0 0 -Y_c*exp(1i*k_c*d) Y_c*exp(-1i*k_c*d);
                   -Y_0*exp(1i*k*L) Y_0*exp(-1i*k*L) Y_0 -Y_0];

               
% encoding the entire j=1 and j=N boundaries

for n = 2:N
    big_ten((1:3) + eqn_count   ,  1:2  ,  (1:4) + (n-2)*2) = unit_ten_j_1;
    eqn_count = eqn_count + 3;
    
    big_ten((1:3) + eqn_count  ,   (2*M-3):(2*M-1)   ,   (1:4) + (n-2)*2) = unit_ten_j_M;
    eqn_count = eqn_count + 3;
end


 

% encoding the known variables
    % pre allocating in-homogeneity vector:
    V = zeros(size(big_ten,1),1);


for j = 1:M
    
    % reflections r from the output:

    big_ten(eqn_count+1 , 2*j-1,2*N) = 1;
    big_ten(eqn_count+1 , 2*j-1,2*N-1) = -1*factor(2);
    eqn_count = eqn_count  +1;
    
    % known t for input / reflections from other inputs
    
    if any(idx_of_input_lines==j) 
        big_ten(eqn_count+1, 2*j-1,1) = 1;
        V(eqn_count+1) = 1;
        
    else
        big_ten(eqn_count+1, 2*j-1,1) = 1;
        big_ten(eqn_count+1 , 2*j-1,2) =  -1*factor(1);
    end
    
    eqn_count = eqn_count  +1;
    
%     setting t_c and r_c for n=1 to zero (in the actual problem they do not exist,
%     and they don't appear in any equation so we can set them to whatever we want)
    if j<M
        big_ten(eqn_count+1 , 2*j,1) = 1;
        eqn_count = eqn_count+1;
        big_ten(eqn_count+1 ,2*j,2) = 1;
        eqn_count = eqn_count+1;
    end
    
   
    
end
toc


%% solving
tic
% translating 3D tensor to 2D matrix:
big_mat = reshape(big_ten, [eqn_count,eqn_count]);

% solving
solution = linsolve(big_mat ,V);
% solution = inv(big_mat)*V;

% getting solution to workable form
solution_mat = reshape(solution, [2*M-1,2*N]);

t = solution_mat(1:2:2*M-1,1:2:2*N-1);
r = solution_mat(1:2:2*M-1,2:2:2*N);


% check energy conservation
cond = abs(sum(abs(t(:,N)).^2) + sum(abs(r(:,1)).^2) - sum(abs(t(:,1)).^2) - sum(abs(r(:,N)).^2) )>1e-12;

if cond
    warning("energy conservation condition does not hold")
end
toc
%% reconstructing voltage and current 

% defining coordinates along the lines:
Npoints = 100; %points per segment
x = linspace(0,N*L,N*Npoints);

% pre-allocating
V = zeros(length(x),M); % voltage
I = zeros(length(x),M); % current


% calculating
for j=1:M
    for n=1:N
       V(Npoints*(n-1)+1:Npoints*n,j) = t(j,n)*exp(1i*k*x(1:Npoints)) + r(j,n)*exp(-1i*k*x(1:Npoints));
       
       I(Npoints*(n-1)+1:Npoints*n,j) = (1/Z_0)*(t(j,n)*exp(1i*k*x(1:Npoints)) - r(j,n)*exp(-1i*k*x(1:Npoints)));
       
    end
end

% calculate power
P = 0.5*real(V.*conj(I));


%% plot - colormap

figure(203)
clf
imagesc(transpose(real(P)), "XData",x )
shading flat
colorbar
yticks(1:M)

title(sprintf("power propagation at %g GHz", Frequency*1e-9))
ylabel( "line" , "fontsize", 15)
xlabel( "position along line (m)" , "fontsize", 15)

colormap jet

%% plot - graphs
figure(204)
clf
plot(x,real(P), 'linewidth', 1.5);
grid on;
xlabel( "position along line (m)" , "fontsize", 15)
ylabel( "power (a.u)" , "fontsize", 15)
title(sprintf("power propagation at %g GHz", Frequency*1e-9),"fontsize", 15)