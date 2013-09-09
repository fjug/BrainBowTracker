function [ graph ] = buildProxGraph( seg, numNuclei, surrpos, surrcol )
%buildProxGraph builds the Proximity Graph from the segmentation data
%   The nodes of the graph are the segmented cells.
%   An edge exists between two nodes iff they are in their surrounding
%   area. Note that the size of this area was defined earlier and this
%   function only takes parameters 'surrpos' and 'surrcol' as inputs in
%   order to determine the graph connections.
%   Note: the graph is stored and returned as a adjacency list...

graph=cell(numNuclei,1); 
msg_n = 0;

fprintf('   buildProxGraph: ');
for nodeFrom = 1:numNuclei
    fprintf(repmat('\b',1,msg_n));
    msg = sprintf('%.1f%%...', 100*(nodeFrom/numNuclei)); % just to be able to monitor progress...
    fprintf('%s', msg);
    msg_n=numel(msg);
    
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
        graph(nodeFrom) = {horzcat(cell2mat(graph(nodeFrom)), nodeTo)}; 
    end
end

fprintf(repmat('\b',1,msg_n));
fprintf('...done!\n');
end

