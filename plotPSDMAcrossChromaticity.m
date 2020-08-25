% function to plot photoreceptor signal distortion metric across chromaticity diagram
% created by ACH 01/07/2020
 
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceSLI.mat');

%% plot for each display (CRT,LCD,Display++)

colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5];
labs = ['S';'M';'L';'R';'I'];

figure('defaultAxesFontSize',18)
for i=1:5
    subplot(3,5,i)
    h(1) = plotDistortionsOverChromaticity(CRT,i);
   
    subplot(3,5,i+5)
    h(2) = plotDistortionsOverChromaticity(LCD,i);
    
    subplot(3,5,10+i)
    h(3) = plotDistortionsOverChromaticity(DP,i);

end
% add in labels
subplot(3,5,6)
ylabel('CIE y');
subplot(3,5,13)
xlabel('CIE x');
subplot(3,5,1)
title('S');
subplot(3,5,2)
title('M');
subplot(3,5,3)
title('L');
subplot(3,5,4)
title('R');
subplot(3,5,5)
title('I');
colorbar;
subplot(3,5,5)

function h = plotDistortionsOverChromaticity(display,d);
    image = reshape(display.meanPhotoreceptorDistortion(d,:,:),[50,50]);
    h=imagesc(image);
    set(gca,'YDir','normal');
    colormap(colorcet('coolwarm'));
    caxis([-60,60]);
    xticks(1:10:50);
    xticklabels(0:0.2:1);
    yticks(1:10:50);
    yticklabels(0:0.2:1);
    axis square

end