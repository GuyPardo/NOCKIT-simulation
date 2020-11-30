
clearvars
v_ph = 1.361104539023962e+06;
v_ph_c = 1.403390237821211e+06;
Z0 = 49.359588669802015;
Z_c = 2.202127166858337e+02;
L_tot_c = 1.569148129658632e-04;
L_tot = 3.626436269560705e-05;




w0 = 1+0.1i;
g = 0.3;

N = 11;
n = N^2;

H = w0*eye(n);

idx = @(i,j) i + N*(j-1);
for i=1:N
    for j = 1:N
        if i<N
        H(idx(i,j), idx(i+1,j)) = g;
        end
        if i>1
        H(idx(i,j), idx(i-1,j)) = g;
        end
        if j<N
        H(idx(i,j), idx(i,j+1)) = g;
        end
        if j>1
        H(idx(i,j), idx(i,j-1)) = g;
        end
    end
end

%%

 w = 0.1:0.01:5;
 T = zeros(length(w),n,n);
 
for i=1:length(w) 
    mat = inv(w(i)*eye(length(H)) - H);
    T(i,:,:) = mat;

    
end

plot(w, abs(T(:,idx(6,1), idx(11,N))).^2)
grid on
    

