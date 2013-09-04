% loading the segmentations and image data at t and (t+1)
fprintf('Starting to load time t...');
if ( exist('seg_t', 'var')~=1 ) 
    vars = load('data/day11/segmentation-1829nuclei-digit2.mat');
    seg_t = vars.ZIBAmira_segmentation_1829nuclei_digit2_mat;
end
if ( exist('img_t_r', 'var')~=1 ) 
    vars = load('data/day11/C1-digit2.tif.mat');
    img_t_r = vars.ZIBAmira_C1_digit2_tif_mat;
end
if ( exist('img_t_g', 'var')~=1 ) 
    vars = load('data/day11/C3-digit2.tif.mat');
    img_t_g = vars.ZIBAmira_C3_digit2_tif_mat;
end
if ( exist('img_t_b', 'var')~=1 ) 
    vars = load('data/day11/C4-digit2.tif.mat');
    img_t_b = vars.ZIBAmira_C4_digit2_tif_mat;
end
fprintf(' ...done!\n');


fprintf('Starting to load time t+1...');
if ( exist('seg_tp1', 'var')~=1 ) 
    vars = load('data/day12/segmentation-day12-1393nuclei-digit2.mat');
    seg_tp1 = vars.ZIBAmira_segmentation_day12_1393nuclei_digit2_mat;
end
if ( exist('img_tp1_r', 'var')~=1 ) 
    vars = load('data/day12/C1-digit2.tif.mat');
    img_tp1_r = vars.ZIBAmira_C1_digit2_tif_mat;
end
if ( exist('img_tp1_g', 'var')~=1 ) 
    vars = load('data/day12/C3-digit2.tif.mat');
    img_tp1_g = vars.ZIBAmira_C3_digit2_tif_mat;
end
if ( exist('img_tp1_b', 'var')~=1 ) 
    vars = load('data/day12/C4-digit2.tif.mat');
    img_tp1_b = vars.ZIBAmira_C4_digit2_tif_mat;
end
fprintf(' ...done!\n');

% setting some metadata
numNuclei_t  =1829;
numNuclei_tp1=1393;
segResX=256;
segResY=512;
maxZ=21;

fprintf('Starting to read segmented cell voxel colors...\n');
if ( exist('colors_t', 'var')~=1 ) 
    colors_t   = getVoxelColors(seg_t,   numNuclei_t,   img_t_r,   img_t_g,   img_t_b);
end
if ( exist('colors_tp1', 'var')~=1 ) 
    colors_tp1 = getVoxelColors(seg_tp1, numNuclei_tp1, img_tp1_r, img_tp1_g, img_tp1_b);
end
fprintf(' ...done!\n');



maxCol=65535;
h=gca;
axis(h,[0 maxCol 0 maxCol 0 maxCol 0 255]);
xlabel(h,'membrane CFP');
ylabel(h,'cytosolic YFP');
zlabel(h,'cytosolic RFP');
dotSize=5;

% mean colors:
figure(1);
format shortG;
for c=1:numNuclei_t
    test=cell2mat(colors_t(c));
    testSize=size(test,2);
    if (testSize>0)
        test=median(double ( cell2mat( colors_t(c) ) ), 2); % double is only needed because MY matlab version sucks... google it...
        testSize=size(test,2);    
        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(test(1,:),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;


figure(2);
format shortG;
% individual color values, or mean+std:
showMeanAndStd=1;
for c=1:numNuclei_t
    test=cell2mat(colors_t(c));
    testSize=size(test,2);
    if (testSize>0)
        if (showMeanAndStd)
            test=mean(cell2mat(colors_t(c)),2);
            dotSize=norm(std(single(cell2mat(colors_t(c))),0,2));
            testSize=size(test,2);    
        end
        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(test(1,:),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;