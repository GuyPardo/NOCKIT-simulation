function [] = plt_network_current(G,coordinates,t_edges,r_edges,freq, dB)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here
res = 60;
if nargin<6
    dB = false;
end


Y = sqrt(G.Edges.C./G.Edges.L);
v_ph = (G.Edges.C.*G.Edges.L).^-0.5;
K = 2*pi*freq./v_ph;

% G.Edges.C.*G.Edges.L
figure(gcf)
hold on

for i = 1:G.numedges
%     l= linspace(0,G.Edges.len(i),res); % coordinate along edge
    x_start = coordinates(1,G.Edges.EndNodes(i,1));
    y_start = coordinates(2,G.Edges.EndNodes(i,1));
    x_end = coordinates(1,G.Edges.EndNodes(i,2));
    y_end = coordinates(2,G.Edges.EndNodes(i,2));
    

    x = linspace(x_start, x_end, res);
    y =linspace(y_start, y_end, res);
    xx = sqrt((x - x_start).^2 + (y-y_start).^2); %% coordinate along line
    z = zeros(size(x));

    current = Y(i)*(t_edges(i)*exp(1i*K(i)*xx)   -  r_edges(i)*exp(-1i*K(i)*xx));
    if dB
        col = 10*log10(abs(real(current)));
    else
        col = real(current);
    end
    
   
    surface([x;x],[y;y],[z;z],[col;col],'facecol','no','edgecol','interp','linew',3);
%     view(2)
   
end



 
 colormap jet;

 colorbar 

end

