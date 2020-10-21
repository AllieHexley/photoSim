% plot Figure 6 of paper
% plots distortions plot note
% i = CRT, ii=LCD, iii=DP
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

% set up colours
CRTcol = [0,0,0];
LCDcol = [0.2,0.2,0.2];
DPcol = [0.6,0.6,0.6];

%% plot fig6ai - CRT distortion Rod

fig = figure('defaultAxesFontSize',12);
plotDistortions(CRT,4,Sim,CRTcol,'A_R','A_R^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6ai.pdf','-dpdf');

%% plot fig6aii - LCD distortion Rod

fig = figure('defaultAxesFontSize',12);
plotDistortions(LCD,4,Sim,LCDcol,'A_R','A_R^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6aii.pdf','-dpdf');

%% plot fig6aii - Display++ distortion Rod

fig = figure('defaultAxesFontSize',12);
plotDistortions(DP,4,Sim,DPcol,'A_R','A_R^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6aiii.pdf','-dpdf');

%% plot fig6bi - CRT distortion Mel

fig = figure('defaultAxesFontSize',12);
plotDistortions(CRT,5,Sim,CRTcol,'A_I','A_I^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6bi.pdf','-dpdf');

%% plot fig6bii - LCD distortion Mel

fig = figure('defaultAxesFontSize',12);
plotDistortions(LCD,5,Sim,LCDcol,'A_I','A_I^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6bii.pdf','-dpdf');

%% plot fig6biii - Display++ distortion Mel

fig = figure('defaultAxesFontSize',12);
plotDistortions(DP,5,Sim,DPcol,'A_I','A_I^*')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6biii.pdf','-dpdf');

%% plot fig6ci - CRT boxplot Rod

fig = figure('defaultAxesFontSize',12);
plotBoxplots(CRT,4,'PSDM_R (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6ci.pdf','-dpdf');

%% plot fig6cii - LCD distortion Rod

fig = figure('defaultAxesFontSize',12);
plotBoxplots(LCD,4,'PSDM_R (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6cii.pdf','-dpdf');

%% plot fig6ciii - DP Rod boxplot

fig = figure('defaultAxesFontSize',12);
plotBoxplots(DP,4,'PSDM_R (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6ciii.pdf','-dpdf');

%% plot fig6di - CRT boxplot mel

fig = figure('defaultAxesFontSize',12);
plotBoxplots(CRT,5,'PSDM_I (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6di.pdf','-dpdf');

%% plot fig6dii - LCD distortion mel

fig = figure('defaultAxesFontSize',12);
plotBoxplots(LCD,5,'PSDM_I (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6dii.pdf','-dpdf');

%% plot fig6diii - DP mel boxplot

fig = figure('defaultAxesFontSize',12);
plotBoxplots(DP,5,'PSDM_I (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,2.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 2];
print(fig, '..\plots\fig6diii.pdf','-dpdf');

%% plot fig6ei - CRT Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(CRT,Sim,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6ei.pdf','-dpdf');

%% plot fig6eii - LCD Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(LCD,Sim,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6eii.pdf','-dpdf');

%% plot fig6eiii - DP Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(DP,Sim,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6eiii.pdf','-dpdf');

%% plot fig6fi - CRT Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(CRT,Sim,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6fi.pdf','-dpdf');

%% plot fig6fii - LCD Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(LCD,Sim,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6fii.pdf','-dpdf');

%% plot fig6fiii - DP Mel

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(DP,Sim,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6fiii.pdf','-dpdf');

%% plot fig6gi - CRT Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(CRT,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6gi.pdf','-dpdf');

%% plot fig6gii - LCD Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(LCD,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6gii.pdf','-dpdf');

%% plot fig6giii - DP Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(DP,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6giii.pdf','-dpdf');

%% plot fig6gi - CRT i

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(CRT,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6hi.pdf','-dpdf');

%% plot fig6hii - LCD i

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(LCD,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6hii.pdf','-dpdf');

%% plot fig6hiii - DP i

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverMB(DP,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig6hiii.pdf','-dpdf');

%%
clear all;

%% functions

% plot spread
function h = plotDistortions(display,d,Sim,col,xlab,ylab)
    h=plot(Sim.ss(d,:),display.ssDistorted(d,:),'.','Color',col);
    hold on;
    plot([0:1],[0:1],'Color',[0.8,0,0.8],'LineWidth',2);
    xlim([0,0.08]);
    ylim([0,0.08]);
    xlabel(xlab);
    ylabel(ylab);
    xticklabels({});
    yticklabels({});
    axis square
    grid on;
    box on;
end

% plot boxplots of spread of distortion
function h = plotBoxplots(display,d,lab)
h=boxplot(display.distortionMetric(d,:),'Orientation','horizontal','Colors','k','Symbol','k.');
set(h,{'linew'},{2});
grid on;
yticklabels({});
xlim([-100,100]);
xlabel(lab);
ylim([0.75,1.25]);
%axis square
box on;
end

% % plot across chromaticity
% function h = plotDistortionsOverChromaticity(display,d,lab);
%     image = reshape(display.meanPhotoreceptorDistortion(d,:,:),[100,100]);
%     h=imagesc(image');
%     set(gca,'YDir','normal');
%     colormap(colorcet('D10'));
%     c=colorbar;
%     c.Label.String = lab;
%     caxis([-60,60]);
%     xticks(1:10:100);
%     xticklabels(0:0.1:1);
%     yticks(1:10:100);
%     yticklabels(0:0.1:1);
%     axis square
%     box on;
%     grid on;
% end

function h = plotDistortionsOverxy(display,Sim,d,lab);
% calculate binned distortions across xy space
% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor signals
% of real-world spectra across chromaticity space

% loop over CIExy diagram split into a 100x100 grid and find mean distortion
% for each photoreceptor of all real-world spectra with a chromaticity that
% falls into a particular grid
xStep=0.01;
yStep=0.01;
i=1; j=1;
x=0:0.01:0.99; y = 0:0.01:0.99;
for i=1:100
    for j=1:100
        xyDistortions(i,j)=mean(display.distortionMetric(d,Sim.xyY(1,:)>x(i) & Sim.xyY(1,:)<x(i)+xStep & Sim.xyY(2,:)>y(j) & Sim.xyY(2,:)<y(j)+yStep));
    end
end
xyDistortions(isnan(xyDistortions)==1)=0;
h=imagesc(xyDistortions');
set(gca,'YDir','normal');
colormap(colorcet('D10'));
c=colorbar;
c.Label.String = lab;
caxis([-60,60]);
xticks(1:10:100);
xticklabels(0:0.1:1);
yticks(1:10:100);
yticklabels(0:0.1:1);
axis square
box on;
grid on;
end

function h=plotDistortionsOverMB(display,d,lab);
% calculate binned distortions across MB space
% inputs:
% 1) structure of display
% 2) d=which photoreceptor distortion to plot
% 3) lab = label of plot

% outputs:
% display structure containing fields for distorted photoreceptor signals
% of real-world spectra across chromaticity space

% loop over CIExy diagram split into a 100x100 grid and find mean distortion
% for each photoreceptor of all real-world spectra with a chromaticity that
% falls into a particular grid
lStep = 0.01;
sStep = 0.01;
mbx = 0:0.01:1;
mby = 0:0.01:1;
for i=1:100
    for j=1:100
        mbDistortions(i,j)=mean(display.distortionMetric(d,display.mbDistorted(2,:)>mbx(i) & display.mbDistorted(2,:)<(mbx(i)+sStep) & display.mbDistorted(1,:)>mby(j) & display.mbDistorted(1,:)<(mby(j)+lStep)),'omitnan'); 
    end
end

mbDistortions(isnan(mbDistortions))=0;
h=imagesc(mbDistortions');
set(gca,'YDir','normal');
colormap(colorcet('D10'));
c = colorbar();
c.Label.String = lab;
c.LineWidth = 2;
zlabel('CIE Y');
caxis([-100,100]);
xticks(1:10:100);
xticklabels(0:0.1:1);
yticks(1:10:100);
yticklabels(0:0.1:1);
ylim([0,20]);
xlim([50,100]);
axis square;
grid on;
end