function [ sub, nodeSet ] = localSubgraph( proxGraph, segNum )
%localSubgraph extract a local subgraph descriptor of the neighborhood of
%cell with segmentation-id 'segNum'.

nodeSet = [];
sub=struct('nodelabels',{[]},'edges',{[]});

% go over all edges and collect all neighboring cell-id's in 'nodeSet'
for edgeId = 1:size(proxGraph.edges,1)
    edge = proxGraph.edges(edgeId,:);
    if ( edge(1)==segNum )
        nodeSet = union( nodeSet, edge(2) );
    end
    if ( edge(2)==segNum )
        nodeSet = union( nodeSet, edge(1) );
    end
end
% create new nodes in subgraph
for nodesetId = 1:size(nodeSet,2)
    sub.nodelabels = vertcat ( sub.nodelabels, proxGraph.nodelabels(nodeSet(nodesetId)) );
end
% and make 'center' node always node 0!
nodeSet = horzcat ( segNum, nodeSet );
sub.nodelabels = vertcat ( proxGraph.nodelabels(segNum), sub.nodelabels );

% go again over edges but now add the ones found before all to 'sub'
for edgeId = 1:size(proxGraph.edges,1)
    edge = proxGraph.edges(edgeId,:);
    nodeFrom = find( nodeSet==edge(1) );
    nodeTo   = find( nodeSet==edge(2) );
    fpTo = edge(3);
    
    if ( edge(1)==segNum           || edge(2)==segNum ) %  || ... 
       % ( ismember(edge(1),nodeSet) && ismember(edge(2),nodeSet) ) )
        if ( isempty(nodeFrom) || isempty(nodeTo) )
            disp('FUCKUP in localSubgraph -- node-id conversion sucky...'); % 'disp' was 'error'
        end
%         fprintf('>>> %4d-->%4d\t\t%4d-->%4d\n',edge(1),edge(2),nodeFrom, nodeTo);
        sub.edges = vertcat(sub.edges, uint32( [nodeFrom nodeTo fpTo] ));
    end
end

end
