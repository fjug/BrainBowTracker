function [] = plotColorsScatter3( colors, labelx, labely, labelz, onlyMeanAndStd )
%plotColorsScatter3 Plots a given <Cx1 cell> containing a list of RGB 
%colors in a 3D scatter plot.
%   colors: a <Cx1 cell> containing RGB voxel colors of segmented objects 
%           as <3xN cell>
%   labelx: x-axis label
%   labely: y-axis label
%   labelz: z-axis label
%   onlyMeanAndStd: if '1' (true), scatter plot contains only one dot per 
%           <3xN cell>, depicting mean color and SD of contained colors.
%           DEFAULT VALUE: 0 (false)

if nargin < 5
    onlyMeanAndStd=0;
end

maxCol=65535;
h=gca;
axis(h,[0 maxCol 0 maxCol 0 maxCol 0 255]);
xlabel(h,labelx);
ylabel(h,labely);
zlabel(h,labelz);
dotSize=5;

format shortG;
% individual color values, or mean+std:
for c=1:size(colors)
    test=cell2mat(colors(c));
    testSize=size(test,2);
    if (testSize>0)
        if (onlyMeanAndStd)
            test=mean(cell2mat(colors(c)),2);
            dotSize=norm(std(single(cell2mat(colors(c))),0,2));
            testSize=size(test,2);    
        end
        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(test(1,:),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;

end

