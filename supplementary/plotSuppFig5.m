% plot Supplementary Figure 5 of paper
% plot PSRM and average PSDM of extra displays
% created by ACH 13/07/2021

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMSRI.mat');

%% and plot the gamuts for these displays
fig = figure('defaultAxesFontSize',12);
plotChromaticity();
hold on;
plot(Sim.xyY(1,:),Sim.xyY(2,:),'w.');
h(1) = plotChromaticityReproduction(Macbook_Pro_2009,[0.5,0.5,0.5]);
h(2) = plotChromaticityReproduction(Macbook_Pro_2014,[0.8,0.8,0.8]);
h(3) = plotChromaticityReproduction(Macbook_Air,[0.2,0.2,0.2]);
h(4) = plotChromaticityReproduction(Surface_Pro,[0.8,0.2,0.2]);
h(5) = plotChromaticityReproduction(NEC,[0.2,0.2,0.8]);
xlabel('CIE x'); ylabel('CIE y');
%legend(h,{'Macbook Pro 2009', 'Macbook Pro 2014','Macbook Air','Surface Pro', 'NEC'});
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, 'supplementary_plots\figS5a.png','-dpng');

%% and bar graph for chromaticity and SMLRI reproduction
fig = figure('defaultAxesFontSize',12);
b = bar([Macbook_Pro_2009.chromaticityReproductionMetric,Macbook_Pro_2009.psdm,Macbook_Pro_2009.realworldReproductionMetric;Macbook_Pro_2014.chromaticityReproductionMetric,Macbook_Pro_2014.psdm,Macbook_Pro_2014.realworldReproductionMetric;Macbook_Air.chromaticityReproductionMetric,Macbook_Air.psdm,Macbook_Air.realworldReproductionMetric;Surface_Pro.chromaticityReproductionMetric,Surface_Pro.psdm,Surface_Pro.realworldReproductionMetric;NEC.chromaticityReproductionMetric,NEC.psdm,NEC.realworldReproductionMetric],'LineWidth',1.5);
b(1).FaceColor='flat'
b(1).CData = [0.3,0,0.6]
b(2).FaceColor='flat'
b(2).CData = [0.6,0,0.3]
b(3).FaceColor='flat'
b(3).CData = [0.3,0.6,0]
xticklabels({'Macbook Pro 2009', 'Macbook Pro 2014','Macbook Air','Surface Pro', 'NEC'});
legend('Chromaticity','PSDM','PSRM');
xlim([0.5,5.5]);
xtickangle(45);
ylabel('Reproduction (%)');
ylim([0,104]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, 'supplementary_plots\figS5b.pdf','-dpdf');

%%
function h = plotChromaticityReproduction(display,col);

hold on;
h=plot(display.xyYMax(1,display.idx),display.xyYMax(2,display.idx),'Color',col,'LineWidth',2);

end




