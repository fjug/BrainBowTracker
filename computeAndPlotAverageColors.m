%day12:
%segmentation=ZIBAmira_segmentation_day12_1393nuclei_digit2_mat;
%membraneSegmentation=ZIBAmira_segmentation_day12_1393nuclei_digit2_mat;
%numNuclei=1393;
%maxZ=24
%day11:
segmentation=ZIBAmira_segmentation_1829nuclei_digit2_mat;
membraneSegmentation=ZIBAmira_segmentation_1829nuclei_digit2_mat;
numNuclei=1829;
maxZ=21;

colors=cell(numNuclei,1); 
membraneIntensity=cell(numNuclei,1); 
for segCoordX= 1:256
    x=(segCoordX-1)*8+4;
    display(x);
for segCoordY=1:512;
    y=(segCoordY-1)*8+4;
for z=1:maxZ
    segCoordZ=(z-1)*4+2;
    whichNucleus=segmentation(segCoordX,segCoordY,segCoordZ);
    if (whichNucleus>0 && whichNucleus<=numNuclei) 
        new = [ZIBAmira_C1_digit2_tif_mat(x,y,z); ZIBAmira_C3_digit2_tif_mat(x,y,z); ZIBAmira_C4_digit2_tif_mat(x,y,z)];
        whichColors = cell2mat(colors(whichNucleus));
        colors(whichNucleus) = {horzcat(whichColors, new)}; 
    end

    whichMembrane = [membraneSegmentation(segCoordX,segCoordY,segCoordZ)];
    if (whichMembrane>0  && whichMembrane<=numNuclei) 
        newMembrane = [ZIBAmira_C1_digit2_tif_mat(x,y,z)];
        whichMembraneIntensity = cell2mat(membraneIntensity(whichMembrane));
        membraneIntensity(whichMembrane) = {horzcat(whichMembraneIntensity, newMembrane)};    
    end
end
end
end


maxCol=65535;
h=gca;
axis(h,[0 maxCol 0 maxCol 0 maxCol 0 255]);
xlabel(h,'membrane CFP');
ylabel(h,'cytosolic YFP');
zlabel(h,'cytosolic RFP');
dotSize=5;

% mean colors, membrane from own image:
figure(1);
format shortG;
for c=1:numNuclei
    test=cell2mat(colors(c));
    testSize=size(test,2);
    if (testSize>0)
        test=median(cell2mat(colors(c)),2);
        testSize=size(test,2);    
        testMembrane=median(cell2mat(membraneIntensity(c)),2);
        str = [c,testMembrane(),test(2,:),test(3,:)];
        display(str);

        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(testMembrane(),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;

figure(3);
format shortG;
% individual color values, or mean+std:
showMeanAndStd=1;
for c=1:numNuclei
    test=cell2mat(colors(c));
    testSize=size(test,2);
    if (testSize>0)
        if (showMeanAndStd)
            test=mean(cell2mat(colors(c)),2);
            dotSize=norm(std(single(cell2mat(colors(c))),0,2));
            %
            testSize=size(test,2);    
        end
        testOnes=ones(1,testSize);
        testCol=testOnes*c;
        scatter3(test(1,:),test(2,:),test(3,:),testOnes*dotSize,testCol);
        hold on;
    end
end
hold off;
