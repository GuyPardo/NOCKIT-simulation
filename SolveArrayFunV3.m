function [t,r] = SolveArrayFunV3(M,N, L,d,k,k_c,Y_0, Y_c,idx_of_input_lines, ref_factor)
% written by Guy 2020_08_26. solves for M lines and N segments per line. 
%   input arguments:
% M - number of lines
% N - number of segments per line (  = one plus the number of couplers per line)
% L - length of each line segment
% d - distance between lines ( = coupling segment length)
% k - wave number for lines
% k_c wave number for couplers
% idx_of_input_lines - an integer or an array of integers in [1,M].
%                      these are the indices of the lines (j) for which we set t(j,1) = 1. 
%  ref_factor - a complex factor that will be scaling the reflections at
%               the boundaries. usually given by (atten^2)*exp(2i*k_coax*coax_l)
%               this should be given as a M*2 matrix: one element for each of the 2M ports
%               the elements corresponding to the port(s) from which the signal is coming will not be used.                
%
%   output arguments:
% t - a MxN matirx containing the t (positively-propagating) amplitudes for the
% main line segments
%  r a MxN matirx containing the r (negatively propagating) amplitudes for the
% main line segments
% this function can be easily modified to also return  t_c and r_c (the
% amplitudes for the couplers). untill now we didn't have a need for that.


%% encoding equations
%  see the pdf for details on the encoding procedure
eqn_num = 4*M*N - 2*M - 2*N +2;


%pre-allocate
big_ten = zeros(eqn_num + 2*M-2, 2*M - 1, 2*N); 
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
    big_ten(eqn_count+1 , 2*j-1,2*N-1) = -1*ref_factor(j,2);
    eqn_count = eqn_count  +1;
    
    % known t for input / reflections from other inputs
    
    if any(idx_of_input_lines==j) 
        big_ten(eqn_count+1, 2*j-1,1) = 1;
        V(eqn_count+1) = 1;
    else
        big_ten(eqn_count+1, 2*j-1,1) = 1;
        big_ten(eqn_count+1 , 2*j-1,2) =   -1*ref_factor(j,1);
    end
    
    eqn_count = eqn_count  +1;
    
%     setting t_c and r_c for n=1 to zero (in the actual problem they do not exist,
%     and they don't appear in any equation so we can set them to whatever we want)
    if j<M
        big_ten(eqn_count+1 , :,:) = full(sparse(2*j,1,[1],size(big_ten,2),size(big_ten,3)));
        eqn_count = eqn_count+1;
        big_ten(eqn_count+1 , :,:) = full(sparse(2*j,2,[1],size(big_ten,2),size(big_ten,3)));
        eqn_count = eqn_count+1;
    end
    
   
    
end


%% solving

% translating 3D tensor to 2D matrix:
big_mat = reshape(big_ten, [eqn_count,eqn_count]);

% solving
solution = linsolve(big_mat ,V);

% getting solution to workable form
solution_mat = reshape(solution, [2*M-1,2*N]);
t = solution_mat(1:2:2*M-1,1:2:2*N-1);
r = solution_mat(1:2:2*M-1,2:2:2*N);


% check energy conservation
cond = abs(sum(abs(t(:,N)).^2) + sum(abs(r(:,1)).^2) - sum(abs(t(:,1)).^2) - sum(abs(r(:,N)).^2) )>1e-12;

if cond
     warning("energy conservation condition does not hold. If there are loss mechanisms then this is to be expected, otherwise it could mean a bug")
end

end

