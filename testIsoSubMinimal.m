G = struct('nodelabels',{[]},'edges',{[]});

G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 65536 0 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 65536 0 0 ] ) );
G.nodelabels = vertcat ( G.nodelabels, encodeFPinUint32( [1 0 0 65536 ] ) );

G.edges = vertcat(G.edges, uint32( [1 2 0 ] )); %G.nodelabels(2)] ));
% G.edges = vertcat(G.edges, uint32( [2 3 0 ] )); %G.nodelabels(2)] ));


perm = [ 2 3 1 ];
[ devnull, inverse_perm ] = sort(perm);
    
[ pG, perm ] = permuteGraphNodes(G, perm);
ppG = permuteGraphNodes(pG, inverse_perm);
g = localSubgraph(G,1);

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
