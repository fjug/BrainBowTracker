function [ A, XY ] = getGraphLayout( G, segpos )
%getGraphLayout Takes VF2 compatible graph data and builds gPlot-input.
%   gPlot is a nice function in order to plot pre-layouted graphs.
%   G is treated as a undirected graph (adj. matrix makes two entries per
%   edge)

n = size(G.nodelabels,1);
m = size(G.edges,1);

A = zeros(n);
XY = zeros(n,2);

if nargin == 2
    for i = 1:n
        positions = segpos{i};
        XY(i,:) = [ positions(1,1) positions(2,1) ];
    end
else
    for i = 1:n
        [ x y ] = pol2cart((2.0*pi*(i-1))/n,1);
        XY(i,1) = x;
        XY(i,2) = y;
    end
end

for i = 1:m
    edge = G.edges(i,:);
    muh = max(edge(3),1);
    A(edge(1),edge(2)) = muh;
    % A(edge(2),edge(1)) = edge(3);
end

end

