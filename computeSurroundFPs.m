function [ sfps, positions, colors ] = computeSurroundFPs( surrDist, xyVoxelSize, zVoxelSize, segpos, seg, segnum, img_r, img_g, img_b )
% computeSurroundFPs Compute finger prints (FPs) in the local surrounding of
% segmented cell areas.
%   surrDist: size of the minimal additional area surrounding the 
%           voxels (in um, as euclidian distance).
%   segpos: list of lists of 3D coordinates giving all voxels corresponding 
%           to i-th cell in the segmentation. I use 'segpos' to compute the
%           lateral extend of the segmented cell and compute the FPs then
%           inside an sphear that is 'radius+latext' wide arround computed
%           center of the segmentation area.
%   seg:    Kainmuellersche segmentation matrix
%   segnum: number of segmented cells (this has to be given because
%           in voxels where two segmentations overlap the cell-numbers
%           are addad up)
%   img_r:  the red channel (x,y,z)
%   img_g:  the green channel (x,y,z)
%   img_b:  the blue channel (x,y,z)

% checking integrity of given params
assert( all(size(img_r)==size(img_g)), 'Given color channel images of unequal size!');
assert( all(size(img_r)==size(img_b)), 'Given color channel images of unequal size!');

% from voxels to real pos (in um)
voxelIndex2RealPos = [xyVoxelSize,xyVoxelSize,zVoxelSize];

% readouts of interest (inits)
sfps     =cell(segnum,1); 
colors   =cell(segnum,1); 
positions=cell(segnum,1); 

% loop over all segmented cells and obtain FP
fprintf('   computeSurroundFPs: ');
msg_n = 0;
for cellId = 1:segnum
    numVoxelsInSurr = 0;
    numFilledVoxelsInSurr = 0;

    fprintf(repmat('\b',1,msg_n));
    msg = sprintf('%.1f%%...', 100*(cellId/segnum)); % just to be able to monitor progress...
    fprintf('%s', msg);
    msg_n=numel(msg);
    
    segVoxels = cell2mat( segpos(cellId) );
    maxCoords = max(segVoxels, [], 2);
    minCoords = min(segVoxels, [], 2);
    voxelSegBoundaries = ceil( (maxCoords-minCoords)/2 );
    centerVoxel = minCoords + voxelSegBoundaries;
    umSegBoundaries = voxelSegBoundaries .* voxelIndex2RealPos';
    umRadius = ceil( norm(umSegBoundaries) ) + surrDist;
    
    % now I simply go over all voxels in a 2*radius sidelength cube
    % centered at centerVoxel and check if distance to centerVoxel<radius
    segXmin = max(1,          round(centerVoxel(1)-umRadius/xyVoxelSize));
    segXmax = min(size(seg,1),round(centerVoxel(1)+umRadius/xyVoxelSize));
    for segX = segXmin:segXmax
        imgCoordX = (segX-1)*8+4;

        segYmin = max(1,          round(centerVoxel(2)-umRadius/xyVoxelSize));
        segYmax = min(size(seg,2),round(centerVoxel(2)+umRadius/xyVoxelSize));
        for segY = segYmin:segYmax
            imgCoordY = (segY-1)*8+4;
            
            segZmin = max(1,          round(centerVoxel(3)-umRadius/zVoxelSize));
            segZmax = min(size(seg,3),round(centerVoxel(3)+umRadius/zVoxelSize));
            for segZ = segZmin:4:segZmax   % remember: seg-z-resolution is currently higher then img-z-res!
                imgCoordZ = 1+floor((segZ-1)/4); % floor? see comment above... ;)
                
                dist = norm([segX-centerVoxel(1),...
                             segY-centerVoxel(2),...
                             segZ-centerVoxel(3)]);
                if (dist < umRadius)
                    numVoxelsInSurr = numVoxelsInSurr + 1;

                    % do not count yourself, but other segmented cell
                    % voxels within radius given (and checked for) above...
                    if (seg(segX,segY,segZ)~=cellId && seg(segX,segY,segZ)>0)
                        % finally we reached the point where we found a voxel
                        % that contributes to the surround FP computation...
                        numFilledVoxelsInSurr = numFilledVoxelsInSurr + 1;

                        newPos = [segX;segY;segZ];
                        whichPositions = cell2mat(positions(cellId));
                        positions(cellId) = {horzcat(whichPositions, newPos)}; 

                        newCol = [img_r(imgCoordX,imgCoordY,imgCoordZ); 
                                  img_g(imgCoordX,imgCoordY,imgCoordZ); 
                                  img_b(imgCoordX,imgCoordY,imgCoordZ)];
                        whichColors = cell2mat(colors(cellId));
                        colors(cellId) = {horzcat(whichColors, newCol)}; 
                    end
                end
            end
        end
    end
                        
    % assemble FP
    voxelCols = cell2mat( colors(cellId) );
    meanCol = mean(voxelCols,2);
    stdCol  = std(double(voxelCols),0,2);
    sfps(cellId) = { vertcat(numVoxelsInSurr,numFilledVoxelsInSurr,meanCol,stdCol) };
end

fprintf(repmat('\b',1,msg_n));
fprintf('...done!\n');
end

