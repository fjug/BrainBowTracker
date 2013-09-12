function [ graph ] = buildProxGraph( seg, numNuclei, surrpos, segFPs )
%buildProxGraph builds the Proximity Graph from the segmentation data where
%all node end edge labels are initialized to 1.
%   The nodes of the graph are the segmented cells.
%   An edge exists between two nodes iff they are in their surrounding
%   area. Note that the size of this area was defined earlier and this
%   function only takes parameters 'surrpos' and 'surrcol' as inputs in
%   order to determine the graph connections.
%   Note: the graph is stored and returned as a adjacency list... layout:
%      g.nodelabels: (n,1) discrete integer labels [L_1 ; L_2 ; ... ; L_n];
%      g.edges: (m,2) edges, [from to] at each line:
%         [e_1_{from} e_1_{to} edgelabel_1 ; ... ; e_m_{from} e_m_{to} edgelabel_m]

graph=struct('nodelabels',{ones(numNuclei,1,'uint32')},'edges',{[]});
msg_n = 0;

fprintf('   buildProxGraph: ');
i=0;
for nodeFrom = 1:numNuclei
    i = i+1;
    
    fprintf(repmat('\b',1,msg_n));
    msg = sprintf('%.1f%%...', 100*(nodeFrom/numNuclei)); % just to be able to monitor progress...
    fprintf('%s', msg);
    msg_n=numel(msg);
    
    % take care of node itself
    graph.nodelabels(i) = encodeFPinUint32( segFPs{i} );
    
    % take care of edges
    positionsToCheck = cell2mat(surrpos(nodeFrom));
    setOfNeighbors = [];
    for surrposIndex = 1:size(positionsToCheck,2)
        nodeTo = seg(positionsToCheck(1,surrposIndex),...
                     positionsToCheck(2,surrposIndex),...
                     positionsToCheck(3,surrposIndex));
        setOfNeighbors = union(setOfNeighbors, nodeTo);
    end
    
    for setIndex = 1:size(setOfNeighbors,2)
        nodeTo = setOfNeighbors(setIndex);
        fpTo = encodeFPinUint32( segFPs{nodeTo} );
        if ( nodeFrom < nodeTo ) 
            graph.edges = vertcat(graph.edges, uint32( [nodeFrom nodeTo fpTo] ));
        end
    end
end

fprintf(repmat('\b',1,msg_n));
fprintf('...done!\n');
end

