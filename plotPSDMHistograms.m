% function to plot hisotgrams showing the spread of the photoreceptor signal distortion metric
% created by ACH 01/07/2020
 
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%% plot for each display (CRT,LCD,Display++)

colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5];
labs = ['S';'M';'L';'R';'I'];

figure('defaultAxesFontSize',18)
for i=1:5
    subplot(3,5,i)
    h(1) = plotHistogram(CRT,i,colsList,1);
   
    subplot(3,5,i+5)
    h(2) = plotHistogram(LCD,i,colsList,2);
    
    subplot(3,5,10+i)
    h(3) = plotHistogram(DP,i,colsList,3);
end

% add in labels
subplot(3,5,6)
ylabel('Frequency');
subplot(3,5,13)
xlabel('Photoreceptor Distortion Metric (%)');    
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
subplot(3,5,5)
l=legend(h,{'CRT','LCD','Display++'});

function h = plotHistogram(display,d,colsList,k)
    h = histogram(display.distortionMetric(d,:),'FaceColor',colsList(k,:)');
    xlim([-50,100]);
    hold on;
    plot([0,0],[0,10000],'m-','LineWidth',2);
    ylim([0,2500]);
    axis square;
end