% function to plot real-world and distorted photoreceptor signals in MacLeod-Boynton space
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
    h(1) = plotDistortions(CRT,i,Sim,colsList,labs,1);
   
    subplot(3,5,i+5)
    h(2) = plotDistortions(LCD,i,Sim,colsList,labs,2);
    
    subplot(3,5,10+i)
    h(3) = plotDistortions(DP,i,Sim,colsList,labs,3);
end
% add in labels
subplot(3,5,6)
ylabel('Distorted Photoreceptor Signal Reproduced');
subplot(3,5,13)
xlabel('Real-World Photoreceptor Signal');
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

function h = plotDistortions(display,d,Sim,colsList,labs,k)
    h=plot(Sim.ss(d,:),display.ssDistorted(d,:),'.','Color',colsList(k,:));
    hold on;
    plot([0:1],[0:1],'b-','LineWidth',2);
    xlim([0,0.12]);
    ylim([0,0.12]);
    title(labs(d));
    xticklabels({});
    yticklabels({});
    axis square
end