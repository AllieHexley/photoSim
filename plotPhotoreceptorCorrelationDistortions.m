% function to plot real-world photoreceptor correlations and photoreceptor correlation distortion metric
% created by ACH 01/07/2020
 
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%% plot for each display (CRT,LCD,Display++)
pairs = [1,2;1,3;1,4;1,5;,2,3;2,4;2,5;3,4;3,5;4,5];
pairNames = ['S','M';'S','L';'S','R';'S','I';'M','L';'M','R';'M','I';'L','R';'L','I';'R','I'];
plotIdx = [2,3,4,5,8,9,10,14,15,20];
colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5];
labs = ['S';'M';'L';'R';'I'];

figure('defaultAxesFontSize',18)
h = plotRealWorldCorrelations(Sim,pairs,plotIdx);

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

function h = plotRealWorldCorrelations(Sim,pairs,plotIdx)
    for i=1:10
        subplot(5,5,plotIdx(i))
        h=plot(Sim.ss(pairs(i,2),:),Sim.ss(pairs(i,1),:),'k.');
        xlim([0,0.12]);
        ylim([0,0.12]);
        axis square
    end
end

function h = plotRealWorldCorrelationHeatmaps(Sim,pairs,plotIdx)
% create matrix of correlations
correlationsTable = [sCorr;mCorr;lCorr;rCorr;iCorr];
% set everything in lower diagonal to 0
i2K = [6,11,12,16,17,18,21,22,23,24];% indexes to keep
for i=1:25
    if ismember(i,i2K)==0
        correlationsTable(i)=0;
    end
end
imagesc(correlationsTable);
hold on;
c=1;
for i=1:5
    for j=1:5
        if ismember(c,i2K)==1
            text(i-0.25,j,['\rho:',num2str(correlationsTable(j,i))],'FontSize',18);
        end
        c=c+1;
    end
end
axis square
myColorMap = parula(256);
%myColorMap(1,:) = 1;
colormap(colorcet('gray'));%myColorMap);
colorbar;
xticks([1:5]);
xticklabels(['S';'M';'L';'R';'I']);
yticks([1:5]);
yticklabels(['S';'M';'L';'R';'I']);
end