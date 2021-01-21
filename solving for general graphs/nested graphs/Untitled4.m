%% define graph
clearvars G
    
L = 100e-6;
d = 20e-6;

M=4;
N=4;


main_nodes = cell(M,N);
secondary_nodes = {};
for i=1:M
    for j=1:N
 
        main_nodes{i,j} = get_id(i,j);
        secondary_nodes = [secondary_nodes, get_id(i,j,'L'),get_id(i,j,'R'),get_id(i,j,'U'),get_id(i,j,'D')];
    end
end
    

G = digraph();
G=G.addnode(main_nodes);
G=G.addnode(secondary_nodes);
% G= G.addedge({'A','B'}, {'F','E'})

%% assign coordinates to nodes

for i=1:M
    for j=1:N
        G.Nodes.X(G.findnode(get_id(i,j))) = (j-1)*L;
        G.Nodes.Y(G.findnode(get_id(i,j))) = -(i-1)*L;
        
        G.Nodes.X(G.findnode(get_id(i,j,'L'))) = (j-1)*L - d;
        G.Nodes.Y(G.findnode(get_id(i,j,'L'))) = -(i-1)*L;
        
        G.Nodes.X(G.findnode(get_id(i,j,'R'))) = (j-1)*L + d;
        G.Nodes.Y(G.findnode(get_id(i,j,'R'))) = -(i-1)*L;
        
        G.Nodes.X(G.findnode(get_id(i,j,'U'))) = (j-1)*L ;
        G.Nodes.Y(G.findnode(get_id(i,j,'U'))) = -(i-1)*L + d;
        
        G.Nodes.X(G.findnode(get_id(i,j,'D'))) = (j-1)*L ;
        G.Nodes.Y(G.findnode(get_id(i,j,'D'))) = -(i-1)*L - d;
        
        
    end
end
%%

% connect each main node to it's secondary nodes
for i=1:M
    for j=1:N
        
        id = get_id(i,j);
        source = {id,id,id,id };
      
        target = {get_id(i,j,'L'),get_id(i,j,'R'),get_id(i,j,'U'),get_id(i,j,'D') };       
        G = G.addedge(source, target, [1,1,1,1]);
        
        
    end
end
clf
plot(G,'xdata', G.Nodes.X, 'ydata', G.Nodes.Y)
%%
for i=2:M-1
    for j=2:N-1
        
        
        source = {get_id(i,j,'L'),get_id(i,j,'R'),get_id(i,j,'U'),get_id(i,j,'D')};
        target = {get_id(i,j-1,'R'),get_id(i,j+1,'L'),get_id(i-1,j,'D'),get_id(i+1,j,'U') };       
        G=G.addedge(source, target, [2,2,2,2]);       

    end
end


% for j=1:N-1
%         source = {get_id(1,j,'L'),get_id(1,j,'R'),get_id(1,j,'U')};
%         target = {get_id(1,j-1,'R'),get_id(1,j+1,'L'),get_id(2,j,'D') };       
%         G=G.addedge(source, target, [2,2,2]);       
% end

figure(1)
clf
plot(G,'xdata', G.Nodes.X, 'ydata', G.Nodes.Y)
grid on
