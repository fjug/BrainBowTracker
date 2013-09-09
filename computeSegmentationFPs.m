function [ fps ] = computeSegmentationFPs(colors)
%computeSegmentationFPs Computes finger prints (FPs) of segmented 
%                       cell interiours.
%   colors: a <Cx1 cell> containing RGB voxel colors of segmented objects 
%           as <3xN cell>

% extract needed sizes here 
% (for better readability below)
segnum = size(colors,1);

% iterate over colors cell array and compute finger print fields
fps=cell(segnum,1);
fprintf('   computeSegmentationFPs: ');
msg_n = 0;
for cellId = 1:segnum
    fprintf(repmat('\b',1,msg_n));
    msg = sprintf('%.1f%%...', 100*(cellId/segnum)); % progress monitor...
    fprintf('%s', msg);
    msg_n=numel(msg);
    
    voxelCols = cell2mat( colors(cellId) );
    colCount = size(voxelCols,2);
    if colCount==0 continue; end % in case color list is empty
    avgC = mean(voxelCols,2);
    stdC = std(double(voxelCols),0,2);
    % assemble FP
    fps(cellId) = { vertcat(colCount,avgC,stdC) };
end

fprintf(repmat('\b',1,msg_n));
fprintf('...done!\n');
end
