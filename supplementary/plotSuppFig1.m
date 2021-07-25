% plot Supplementary Figure 1 of paper
% plot display primaries for MPHDR and Manchester display
% created by ACH 13/07/2021

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('supplementary_data/MPHDR.mat');
manchester = readtable('supplementary_data/Manchester_5Primary.xlsx');
manchester_primaries = table2array(manchester(:,2:6));
manchester_wls = table2array(manchester(:,1));
load('photosimMetrics_ReproduceLMSRI.mat');

%% plot primaries for Manchester display

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(manchester_wls,manchester_primaries(:,5),'Color',[0.5,0,0],'LineWidth',2);
h(2)=plot(manchester_wls,manchester_primaries(:,3),'Color',[0,0.5,0],'LineWidth',2);
h(3)=plot(manchester_wls,manchester_primaries(:,2),'Color',[0,0,0.5],'LineWidth',2);
h(3)=plot(manchester_wls,manchester_primaries(:,1),'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h(3)=plot(manchester_wls,manchester_primaries(:,4),'Color',[0.9290,0.6940,0.1250],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,700]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, 'supplementary_plots\figS1a.pdf','-dpdf');

%% plot primaries for MPHDR

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(wls,mphdrRGBCMY(:,1),'Color',[0.5,0,0],'LineWidth',2);
h(2)=plot(wls,mphdrRGBCMY(:,4),'Color',[0.8500, 0.3250, 0.0980],'LineWidth',2);
h(3)=plot(wls,mphdrRGBCMY(:,6),'Color',[0,0,0.5],'LineWidth',2);
h(3)=plot(wls,mphdrRGBCMY(:,3),'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h(3)=plot(wls,mphdrRGBCMY(:,2),'Color',[0.9290,0.6940,0.1250],'LineWidth',2);
h(3)=plot(wls,mphdrRGBCMY(:,5),'Color',[0,0.5,0],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,700]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, 'supplementary_plots\figS1b.pdf','-dpdf');

%% and plot the gamuts for these displays
fig = figure('defaultAxesFontSize',12);
plotChromaticity();
hold on;
plot(Sim.xyY(1,:),Sim.xyY(2,:),'w.');
h(1) = plotChromaticityReproduction(Man,[0.5,0.5,0.5]);
h(2) = plotChromaticityReproduction(MPHDR,[0.8,0.8,0.8]);
xlabel('CIE x'); ylabel('CIE y');
legend(h,{'Manchester VDU', 'Oxford MPHDR'});
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, 'supplementary_plots\figS1c.pdf','-dpdf');

%% and bar graph for chromaticity and SMLRI reproduction
fig = figure('defaultAxesFontSize',12);
b = bar([Man.chromaticityReproductionMetric,Man.realworldReproductionMetric;MPHDR.chromaticityReproductionMetric,MPHDR.realworldReproductionMetric],'LineWidth',1.5);
b(1).FaceColor='flat'
b(1).CData = [0.3,0,0.6]
b(2).FaceColor='flat'
b(2).CData = [0.6,0,0.3]
xticklabels({'5P-VDU','MPHDR'});
%legend('Chromaticity','PSRM');
xlim([0.5,2.5]);
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
print(fig, 'supplementary_plots\figS1d.pdf','-dpdf');

%%
function h = plotChromaticityReproduction(display,col);

hold on;
h=plot(display.xyYMax(1,display.idx),display.xyYMax(2,display.idx),'Color',col,'LineWidth',2);

end




