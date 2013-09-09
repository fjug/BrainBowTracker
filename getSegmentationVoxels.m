function [ positions, colors ] = getSegmentationVoxels( seg, segnum, img_r, img_g, img_b )
% reading positions and colors of all voxels of segmented cells from image data
%   seg:    Kainm?llersche segmentation matrix
%   segnum: number of segmented cells (this has to be given because
%           in voxels where two segmentations overlap the cell-numbers
%           are addad up)
%   img_r:  the red channel (x,y,z)
%   img_g:  the green channel (x,y,z)
%   img_b:  the blue channel (x,y,z)

% checking integrity of given params
assert( all(size(img_r)==size(img_g)), 'Given color channel images of unequal size!');
assert( all(size(img_r)==size(img_b)), 'Given color channel images of unequal size!');

% extract needed sizes here 
% (for better readability below)
segResX = size(seg,1);
segResY = size(seg,2);
maxZ = size(img_r,3);

% iterate over segmentation matrix and look up values in color channels
% % GoodToKnows:
% (o) resolution of segmentation might be different from image resolution
% (o) currently: x/y resolutions in segmentations reduced
% (o) z-resolution augmented
colors=cell(segnum,1); 
positions=cell(segnum,1); 
msg_n = 0;

fprintf('   getSegmentedVoxels: ');
for segCoordX = 1:segResX
    fprintf(repmat('\b',1,msg_n));
    msg = sprintf('%.1f%%...', 100*(segCoordX/segResX)); % just to be able to monitor progress...
    fprintf('%s',msg);
    msg_n=numel(msg);
    
    imgCoordX = (segCoordX-1)*8+4;
    
    for segCoordY = 1:segResY
        imgCoordY = (segCoordY-1)*8+4;
        for imgCoordZ = 1:maxZ
            segCoordZ = (imgCoordZ-1)*4+2;
            whichNucleus=seg(segCoordX,segCoordY,segCoordZ);
            if (whichNucleus>0 && whichNucleus<=segnum) 
                newPos = [segCoordX;segCoordY;segCoordZ];
                whichPositions = cell2mat(positions(whichNucleus));
                positions(whichNucleus) = {horzcat(whichPositions, newPos)}; 
        

                newCol = [img_r(imgCoordX,imgCoordY,imgCoordZ); 
                          img_g(imgCoordX,imgCoordY,imgCoordZ); 
                          img_b(imgCoordX,imgCoordY,imgCoordZ)];
                whichColors = cell2mat(colors(whichNucleus));
                colors(whichNucleus) = {horzcat(whichColors, newCol)}; 
            end
        end
    end
end

fprintf(repmat('\b',1,msg_n));
fprintf('...done!\n');
end

