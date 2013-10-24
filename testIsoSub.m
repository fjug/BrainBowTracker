G = struct('nodelabels',{[]},'edges',{[]});

% 1
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 65536 0 ] ) );
% 2-5
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 65536 0 0 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 0 65536 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 65536/2 65536/2 65536/2 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 65536 0 ] ) );
%6-8
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 0 0 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 0 0 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 0 0 ] ) );

G.edges = vertcat(G.edges, uint32( [1 2 0 ] )); %G.nodelabels(2)] ));
G.edges = vertcat(G.edges, uint32( [1 3 0 ] )); %G.nodelabels(3)] ));
G.edges = vertcat(G.edges, uint32( [5 2 0 ] )); %G.nodelabels(2)] ));
G.edges = vertcat(G.edges, uint32( [5 3 0 ] )); %G.nodelabels(3)] ));
G.edges = vertcat(G.edges, uint32( [1 6 0 ] )); %G.nodelabels(6)] ));
G.edges = vertcat(G.edges, uint32( [2 7 0 ] )); %G.nodelabels(7)] ));
G.edges = vertcat(G.edges, uint32( [5 4 0 ] )); %G.nodelabels(4)] ));
G.edges = vertcat(G.edges, uint32( [3 8 0 ] )); %G.nodelabels(8)] ));

perm = [ 6     8     5     4     2     1     7     3 ];
[ devnull, inverse_perm ] = sort(perm);
    
[ pG, perm ] = permuteGraphNodes(G, perm);
ppG = permuteGraphNodes(pG, inverse_perm);
g = localSubgraph(G,5);

figure(1);
[ gA, gXY ] = getGraphLayout(g);
gplot(gA,gXY,'-*');

figure(2);
[ GA, GXY ] = getGraphLayout(G);
gplot(GA,GXY,'-*');

figure(3);
[ pGA, pGXY ] = getGraphLayout(pG);
gplot(pGA,pGXY,'-*');

figure(4);
[ ppGA, ppGXY ] = getGraphLayout(ppG);
gplot(ppGA,ppGXY,'-*');

[count,mappings]   = graphmatch (g, G, 0, 0); count
[count2,mappings2] = graphmatch (g, pG, 0, 0); count2
[count3,mappings3] = graphmatch (g, ppG, 0, 0); count3
