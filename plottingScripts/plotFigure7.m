% plot Figure 7 of paper
% plots correlation distortions plot
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% plot fig7a - CRT correlation distortions

fig = figure('defaultAxesFontSize',12);
h = plotCorrelationDistortionsMatrix(CRT, 'PCDM (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig7a.pdf','-dpdf');

%% plot fig7a - LCD correlation distortions

fig = figure('defaultAxesFontSize',12);
h = plotCorrelationDistortionsMatrix(LCD, 'PCDM (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig7b.pdf','-dpdf');

%% plot fig7a - DP correlation distortions

fig = figure('defaultAxesFontSize',12);
h = plotCorrelationDistortionsMatrix(DP, 'PCDM (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig7c.pdf','-dpdf');

%%
clear all;

%% functions
function h = plotCorrelationDistortionsMatrix(display,lab)
correlationsTable = [0,display.photoreceptorCorrelationsDistortion(1),display.photoreceptorCorrelationsDistortion(2),display.photoreceptorCorrelationsDistortion(3),display.photoreceptorCorrelationsDistortion(4);...
    0, 0, display.photoreceptorCorrelationsDistortion(5), display.photoreceptorCorrelationsDistortion(6), display.photoreceptorCorrelationsDistortion(7); ...
    0, 0, 0, display.photoreceptorCorrelationsDistortion(8), display.photoreceptorCorrelationsDistortion(9);...
    0, 0, 0, 0, display.photoreceptorCorrelationsDistortion(10)];
i2K = [6,11,12,16,17,18,21,22,23,24];% indexes to keep
h=imagesc((correlationsTable));
hold on;
c=1;
for i=1:5
    for j=1:5
        if ismember(c,i2K)==1
            text(i-0.3,j,[num2str(round(correlationsTable(j,i),3))],'FontSize',8);
        end
        c=c+1;
    end
end
colormap(colorcet('D10'));
%caxis([-20,20]);
caxis([-0.1,0.1]);
c = colorbar;
c.Label.FontSize = 12;
c.Label.String = lab;
c.LineWidth = 2;
axis square
set(gca,'XAxisLocation','top','YAxisLocation','right')
xticks([1:5]);
xticklabels(['S';'M';'L';'R';'I']);
yticks([1:5]);
yticklabels(['S';'M';'L';'R';'I']);
end