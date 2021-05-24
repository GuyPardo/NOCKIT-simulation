function [] = plt_nockit_current(G,coordinates,t_edges,r_edges,freq, dB)
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

ax1 = axes;
hold on
ax2 = axes;
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
    
   if G.Edges.Weight(i)==1
    surface(ax1,[x;x],[y;y],[col;col],'facecol','no','edgecol','interp','linew',4);

   end
   if G.Edges.Weight(i)==2
    surface(ax2,[x;x],[y;y],[col;col],'facecol','no','edgecol','interp','linew',8);
   end
   
   
end


linkaxes([ax1,ax2])
ax2.Visible = 'off';
ax2.XTick = [];
ax2.YTick = [];
ax1.Visible = 'off';
ax1.XTick = [];
ax1.YTick = [];

colormap(ax1,'copper');
colormap(ax2, 'winter');
% % set([ax1,ax2],'Position',[.17 .11 .685 .815]);
cb1 = colorbar(ax1,'Position',[.09 .11 .0475 .815]  );
cb2 = colorbar(ax2,'Position',[.82 .11 .0475 .815]);
cb1.Label.String = 'coupler current (A)';
cb1.FontSize = 14;
cb2.Label.String = 'traces current (A)';
cb2.FontSize = 14;

%  colormap jet;
% 
%  colorbar 


end

