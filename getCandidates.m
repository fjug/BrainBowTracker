function [ candidates ] = getCandidates( mappings )
% getCandidates 
%   Extracts a a set of numbers from the first column of given matrix

candidates = [];

for i=1:size(mappings,1)
    candidates = union ( candidates, mappings(i,1) );
end