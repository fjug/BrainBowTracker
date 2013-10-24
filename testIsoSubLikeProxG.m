delete('scheiss.log');

G = struct('nodelabels',{[]},'edges',{[]});

% 1
G.nodelabels = vertcat ( G.nodelabels, uint32(15) );
% 2-5
G.nodelabels = vertcat ( G.nodelabels, uint32(58) );
G.nodelabels = vertcat ( G.nodelabels, uint32(79) );
G.nodelabels = vertcat ( G.nodelabels, uint32(183) );
G.nodelabels = vertcat ( G.nodelabels, uint32(191) );

%central star
G.edges = vertcat(G.edges, uint32( [1 2 0 ] )); 
G.edges = vertcat(G.edges, uint32( [1 3 0 ] )); 
G.edges = vertcat(G.edges, uint32( [1 4 0 ] )); 
G.edges = vertcat(G.edges, uint32( [1 5 0 ] )); 
G.edges = vertcat(G.edges, uint32( [5 1 0 ] )); 
G.edges = vertcat(G.edges, uint32( [4 1 0 ] )); 
G.edges = vertcat(G.edges, uint32( [3 1 0 ] )); 
G.edges = vertcat(G.edges, uint32( [2 1 0 ] )); 

%additional inner edges
G.edges = vertcat(G.edges, uint32( [2 4 0 ] )); 
G.edges = vertcat(G.edges, uint32( [4 2 0 ] )); 



g = localSubgraph(G,1);

figure(1);
[ gA, gXY ] = getGraphLayout(g);
gplot(gA,gXY,'-*');

figure(2);
[ GA, GXY ] = getGraphLayout(G);
gplot(GA,GXY,'-*');

% figure(3);
% [ gA, gXY ] = getGraphLayout(g);
% gplot(gA,gXY,'-*');
% 
% figure(4);
% [ GA, GXY ] = getGraphLayout(G);
% gplot(GA,GXY,'-*');

[count,mappings]   = graphmatch (g, G, 1, 1); count
