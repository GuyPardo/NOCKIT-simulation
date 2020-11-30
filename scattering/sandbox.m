% simulating response of two coupled resonators with resnanat frequency w0
% and coupling g
clearvars
Nmax = 1;

a = spdiags(sqrt((1:Nmax)'),-1,Nmax+1,Nmax+1)  ;%anihilation 
a_dag = a';%creation 


w1 = 1+0i;
w2 = 2+0i;
g =0.5 ;


 H = [w1,g;g,w2]; 
 
 
 %% calculate frequency response with the lippman-schwinger-ish formula
 
 w = 0.1:0.01:20;
 T1 = zeros(size(w));
 T2 = zeros(size(w));
for i=1:length(w) 
    mat = inv(w(i)*eye(length(H)) - H);
    T1(i) = mat(1,1);
    T2(i) = mat(2,2);
    
end
    
 plot(w,abs(T1).^2)
 
 %%
 eigs(real(H))
 