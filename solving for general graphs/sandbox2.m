tic
clearvars
v_ph = 1.3611e+06;
v_ph_c =  1.4088e+06;
Y0 = 0.0203;
Yc = 0.0027;

N=30;
freq = 6e9;
L0 = 100e-6;
d = 20e-6;
k0 = 2*pi*freq/v_ph;
kc = 2*pi*freq/v_ph_c;




e = ones(2*N+4,1);


% A = spdiags([ e 2*e], [1,2], 2*N+4, 2*N+4);
% A(1,2)=0; A(2,1) = 0; A(2*N+3, 2*N+4) =0; A(2*N+4, 2*N+3) =0;
% A(2,3) = 0; A(3,2) =0; A(2*N+2, 2*N+3) =0; A(2*N+3, 2*N+2) =0;

% define ladder grpah:
source = [3:2:2*N+1, 1:2:2*N+1, 2:2:2*N+2];
target = [4:2:2*N+2, 3:2:2*N+3, 4:2:2*N+4];
w = [1*ones(1,N), 2*ones(1,length(source)-N)]; % 1 is for coupling edges, 2 for main edges

G = digraph(source,target,w);
edge_num = G.numedges;
%clearvars G.Edges.k G.Edges.L G.Edges.Y
G.Edges.k(G.Edges.Weight==2) = k0*ones(sum(G.Edges.Weight==2),1);
G.Edges.k(G.Edges.Weight==1) = kc*ones(sum(G.Edges.Weight==1),1);
G.Edges.L(G.Edges.Weight==2) = L0*ones(sum(G.Edges.Weight==2),1);
G.Edges.L(G.Edges.Weight==1) = d*ones(sum(G.Edges.Weight==1),1);
G.Edges.Y(G.Edges.Weight==2) = Y0*ones(sum(G.Edges.Weight==2),1);
G.Edges.Y(G.Edges.Weight==1) = Yc*ones(sum(G.Edges.Weight==1),1);

disp('graph construction:')
toc

figure(504)
G.plot('EdgeLabel', G.Edges.Weight);

tIdx = @(x) 2*x-1;
rIdx = @(x) 2*x;
%% encode equations
tic
mat = zeros(2*edge_num);
vec = zeros(2*edge_num,1);
eqn_count = 0;


% loop on nodes
for i=1:G.numnodes
    outedges = G.outedges(i);
    inedges = G.inedges(i);
    
    edges = [ inedges; outedges];
    
    % edges parameters
    L_in = G.Edges.L(inedges);
    L_out = zeros(size(outedges));
    
    k_out = G.Edges.k(outedges);
    Y_out  = G.Edges.Y(outedges);
    
    k_in = G.Edges.k(inedges);
    Y_in = G.Edges.Y(inedges);
    
    L = [L_in; L_out];
    k = [k_in; k_out];
    Y = [Y_in; -Y_out];  % the minus sign is to differentiate btw currnt into the node and out of the nod
    
    
    % encode voltage continuity equation for the  ith node: loop on node edges
    for j = 1:length(edges)-1
       edge = edges(j);
       next_edge = edges(j+1);
        % voltage eqn:
        mat(eqn_count+1,tIdx(edge)) = exp(1i*L(j)*k(j)); 
        mat(eqn_count+1,rIdx(edge)) = exp(-1i*L(j)*k(j));
        mat(eqn_count+1,tIdx(next_edge)) = -exp(1i*L(j+1)*k(j+1));
        mat(eqn_count+1,rIdx(next_edge)) = -exp(-1i*L(j+1)*k(j+1));
        eqn_count = eqn_count+1; 

    end
    
    % encode the ith node current equation:
    if length(edges)>1
        mat(eqn_count+1, tIdx(edges)) = Y.*exp(1i*k.*L);
        mat(eqn_count+1, rIdx(edges)) = -Y.*exp(-1i*k.*L);
        eqn_count = eqn_count+1;
    end
 
end

% boundary conditions:

mat(eqn_count+1, tIdx(G.findedge(1,3)))=1;
vec(eqn_count+1) = 1;
eqn_count = eqn_count+1;

mat(eqn_count+1, tIdx(G.findedge(2,4)))=1; eqn_count = eqn_count+1;
mat(eqn_count+1, rIdx(G.findedge(2*N+2,2*N+4)))=1; eqn_count = eqn_count+1;
mat(eqn_count+1, rIdx(G.findedge(2*N+1,2*N+3)))=1; eqn_count = eqn_count+1;
disp('encode eqns:')
toc
%% solve
tic
solution = linsolve(mat,vec);

t_edges = solution(1:2:end-1);
r_edges = solution(2:2:end);

 t = [transpose(t_edges(G.findedge(1:2:2*N+1,3:2:2*N+3))); transpose(t_edges(G.findedge(2:2:2*N+2,4:2:2*N+4)))];
 r = [transpose(r_edges(G.findedge(1:2:2*N+1,3:2:2*N+3))); transpose(r_edges(G.findedge(2:2:2*N+2,4:2:2*N+4)))];
 
 % check energy conservation:
 abs(r(1,1))^2 + abs(r(2,1))^2+ abs(t(1,end))^2 + abs(t(2,end))^2
disp('solve:')
toc
 %% reconstruct physical quantities
% defining coordinates along the lines:
Npoints = 100; %points per segment
M=2;
x = linspace(0,(N+1)*L0,(N+1)*Npoints);

% pre-allocating
V = zeros(length(x),M); % voltage
I = zeros(length(x),M); % current


% calculating
for j=1:M
    for n=1:N+1
       V(Npoints*(n-1)+1:Npoints*n,j) = t(j,n)*exp(1i*k0*x(1:Npoints)) + r(j,n)*exp(-1i*k0*x(1:Npoints));
       
       I(Npoints*(n-1)+1:Npoints*n,j) = Y0*(t(j,n)*exp(1i*k0*x(1:Npoints)) - r(j,n)*exp(-1i*k0*x(1:Npoints)));
       
    end
end

% calculate power
P = 0.5*real(V.*conj(I));


%% plot
figure(502)
plot(x,P, 'linewidth', 1.5); grid on; legend('line 1', 'line 2','fontsize',14, 'location', 'best');

xlabel( "position along line (m)" , "fontsize", 15)
ylabel( "power (a.u.)" , "fontsize", 15)
title(sprintf("power propagation at %g GHz", freq*1e-9),"fontsize", 15)