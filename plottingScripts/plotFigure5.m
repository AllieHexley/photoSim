% plot Figure 5 of paper
% plots real-world and chromaticity reproduction plots
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%% plot fig5a - Chromaticity diagram % capture

fig = figure('defaultAxesFontSize',12);
plotChromaticity();
hold on;
plot(Sim.xyY(1,:),Sim.xyY(2,:),'w.');
h(1)=plotChromaticityReproduction(CRT,[0,0,0]);
h(2)=plotChromaticityReproduction(LCD,[0.2,0.2,0.2]);
h(3)=plotChromaticityReproduction(DP,[0.6,0.6,0.6]);
h(4)=plotChromaticityReproduction(FP1,[0,0.8,0]);
h(5)=plotChromaticityReproduction(FP2,[0,0,0.8]);
xlabel('CIE x'); ylabel('CIE y');
%legend(h,{'CRT','LCD','Display++','Narrowband 5P', 'Broadband 5p'});
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig5a.pdf','-dpdf');

%% plot fig5b - Reproduction of chromaticity bar graph

fig = figure('defaultAxesFontSize',12);
b = bar([CRT.chromaticityReproductionMetric,LCD.chromaticityReproductionMetric,DP.chromaticityReproductionMetric,FP1.chromaticityReproductionMetric,FP2.chromaticityReproductionMetric],'LineWidth',1.5,'FaceColor',[0.3,0,0.6]);
xticklabels({'CRT','Dell LCD','Display++', 'NB 5P', 'BB 5P'});
xlim([0.5,5.5]);
xtickangle(45);
ylabel('Chromaticity Reproduction (%)');
ylim([0,100]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig5b.pdf','-dpdf');

%% plot fig5c - Reproduction of full photoreceptor signals bar graph

fig = figure('defaultAxesFontSize',12);
b = bar([CRT.realworldReproductionMetric,LCD.realworldReproductionMetric,DP.realworldReproductionMetric,FP1.realworldReproductionMetric,FP2.realworldReproductionMetric],'LineWidth',1.5,'FaceColor',[0.3,0,0.6]);
xticklabels({'CRT','Dell LCD','Display++', 'NB 5P', 'BB 5P'});
xlim([0.5,5.5]);
xtickangle(45);
ylabel('PSRM (%)');
ylim([0,100]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig5c.pdf','-dpdf');

%% functions

function h = plotChromaticityReproduction(display,col);

hold on;
h=plot(display.xyYMax(1,display.idx),display.xyYMax(2,display.idx),'Color',col,'LineWidth',2);

end
