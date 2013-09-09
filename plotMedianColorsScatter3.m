function [] = plotMedianColorsScatter3( colors, labelx, labely, labelz )
%plotMedianColorScatter3 Plots a given <Cx1 cell> containing a list of RGB 
%colors in a 3D scatter plot.
%   colors: a <Cx1 cell> containing RGB voxel colors of segmented objects 
%           as <3xN cell>
%   labelx: x-axis label
%   labely: y-axis label
%   labelz: z-axis label

maxCol=65535;
h=gca;
axis(h,[0 maxCol 0 maxCol 0 maxCol 0 255]);
xlabel(h,labelx);
ylabel(h,labely);
zlabel(h,labelz);
dotSize=5;

% mean colors:
format shortG;
for c=1:size(colors)
    test=cell2mat(colors(c));
    testSize=size(test,2);
    if (testSize>0)
        test=median(double ( cell2mat( colors(c) ) ), 2); % double is only needed because MY matlab version sucks... google it...
        testSize=size(test,2);    
        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(test(1,:),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;

end

