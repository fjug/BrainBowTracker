function [ newG, perm ] = permuteGraphNodes( oldG, perm )
%permuteGraphNodes Guess from title... what might happen here?

newG = struct('nodelabels',{[]},'edges',{[]});
if nargin < 2
    perm = randperm ( size(oldG.nodelabels, 1) );
end
[ devnull, inverse_perm ] = sort(perm);

newG.nodelabels = oldG.nodelabels(perm);

for edgeId = 1:size(oldG.edges,1)
    edge = oldG.edges(edgeId,:);
    newG.edges = vertcat(newG.edges, uint32( [inverse_perm(edge(1))...
                                              inverse_perm(edge(2))...
                                              edge(3)] ));
end

end

