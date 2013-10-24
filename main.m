addpath('./gboost-0.1.1/bin/');
delete('scheiss.log');

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
numNuclei_t   = 1829;
numNuclei_tp1 = 1393;
segResX=256;
segResY=512;
maxZ=21;

fprintf('Starting to read segmented cell voxel colors...');
do = 1;
if ( exist('segpos_t', 'var')~=1 || exist('colors_t', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    [segpos_t,   colors_t]   = getSegmentationVoxels(seg_t,   numNuclei_t,   img_t_r,   img_t_g,   img_t_b);
end
if ( exist('segpos_tp1', 'var')~=1 || exist('colors_tp1', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    [segpos_tp1, colors_tp1] = getSegmentationVoxels(seg_tp1, numNuclei_tp1, img_tp1_r, img_tp1_g, img_tp1_b);
end
fprintf(' ...done!\n');


if (0)
    % median color per cell
    figure(1);
    fprintf('Plotting median colors for t...');
    plotMedianColorsScatter3(colors_t,'membrane CFP','cytosolic YFP','cytosolic RFP');
    fprintf(' ...done!\n');

    % mean and sd of colors per cell
    figure(2);
    fprintf('Plotting colors for t...');
    plotColorsScatter3(colors_t,'membrane CFP','cytosolic YFP','cytosolic RFP',1);
    fprintf(' ...done!\n');

    % all individual voxel colors
    figure(3);
    fprintf('Plotting colors for t...');
    plotColorsScatter3(colors_t,'membrane CFP','cytosolic YFP','cytosolic RFP',0);
    fprintf(' ...done!\n');
end

fprintf('Starting fingerprinting...');
do = 1;
% compute finger prints (FP's) of segmented cell areas
if ( exist('segFPs_t', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    segFPs_t   = computeSegmentationFPs(colors_t);
end
if ( exist('segFPs_tp1', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    segFPs_tp1 = computeSegmentationFPs(colors_tp1);
end
% compute finger prints (FP's) of segmented cell surrounds
if ( exist('surrFPs_t', 'var')~=1 || exist('surrpos_t', 'var')~=1 || exist('surrcol_t', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    [surrFPs_t,surrpos_t,surrcol_t] = computeSurroundFPs(10.0, 0.2768*8, 10.0/4,... 
                                                         segpos_t, seg_t, numNuclei_t,...
                                                         img_t_r, img_t_g, img_t_b);
end
if ( exist('surrFPs_tp1', 'var')~=1 || exist('surrpos_tp1', 'var')~=1 || exist('surrcol_tp1', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    [surrFPs_tp1,surrpos_tp1,surrcol_tp1] = computeSurroundFPs(20.0, 0.2768*8, 10.0/4,... 
                                                         segpos_tp1, seg_tp1, numNuclei_tp1,...
                                                         img_tp1_r, img_tp1_g, img_tp1_b);
end
fprintf(' ...done!\n');

fprintf('Starting building ProxGraph...');
do = 1;
% build a graph out of the segmentation 
% (cells are connected if they are close to each other)
if ( exist('proxGraph_t', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    proxGraph_t = buildProxGraph(seg_t, numNuclei_t, surrpos_t, segFPs_t);
end
if ( exist('proxGraph_tp1', 'var')~=1 ) 
    if do
        fprintf('\n');
        do=0;
    end
    proxGraph_tp1 = buildProxGraph(seg_tp1, numNuclei_tp1, surrpos_tp1, segFPs_tp1);
end
fprintf(' ...done!\n');

figure(1);
[ GA, GXY ] = getGraphLayout(proxGraph_t,segpos_t);
gplot(GA,GXY,'-*');
figure(2);
[ GA, GXY ] = getGraphLayout(proxGraph_tp1,segpos_tp1);
gplot(GA,GXY,'-*r');


min_col_slack = 3;
max_col_slack = 32;
min_subg_n = 4;
max_mappings = 1;
candidates = [];
for center=[ 8 38 162 1483 781 ]
    fprintf('\nCELL: %d\n==========\n', center);
    
    % cut local subgraph centered at current cell to look at
    subg = localSubgraph( proxGraph_t, center );

    % and show the beast!
    figure(3);
    [ gA, gXY ] = getGraphLayout(subg,segpos_t);
    gplot(gA,gXY,'-*');
    figure(4);
    [ gA, gXY ] = getGraphLayout(subg);
    gplot(gA,gXY,'-*');
    
    % check if subg is large enough to start search in t+1
    if ( length(subg.nodelabels) < min_subg_n ) 
        fprintf('Local subgraph is below minimal size for search to make sense...\n');
        continue;
    end
    
    % now start looking for monomorphisms
    % note: color_slack is set to min at first and gets increased until a
    % match is found. This could/should be replaced by some binary search!
    color_slack = min_col_slack;
    % THE MATCHING
    [count,mappings] = graphmatch (subg, proxGraph_tp1, max_mappings, 1, 255, color_slack);
    
    while ( count == 0 && color_slack <= max_col_slack )
        fprintf('No results using a color-slack of %d...\n',color_slack);
        color_slack = color_slack + 1;
        [count,mappings] = graphmatch (subg, proxGraph_tp1, max_mappings, 1, 255, color_slack);
    end

    if ( count > 0 )
        fprintf('Cell %d -- %d of at max %d MATCH(ES) FOUND!!!\n', center, count, max_mappings);
        candidates = getCandidates(mappings);
        candidates
        
        c1 = candidates(1);
        matchedG = localSubgraph( proxGraph_tp1, c1 );
        
        figure(5);
        [ gA, gXY ] = getGraphLayout(matchedG,segpos_tp1);
        gplot(gA,gXY,'-*r');
        figure(6);
        [ gA, gXY ] = getGraphLayout(matchedG);
        gplot(gA,gXY,'-*r');
    else
        fprintf('Cell %d -- Negativo!\n',center);
    end
    
    fprintf('\npress any key to continue...\n');
    pause;
end