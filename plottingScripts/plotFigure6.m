% plot Figure 6
% plots distortions plot note
% i = CRT, ii=LCD, iii=DP
% created by ACH 01/07/2020

%% load data
%clear all;
close all;
clc;
     
%% load relevant data file

load('../photosimMetrics_ReproduceLMS.mat');

% set up colours
CRTcol = [0,0,0];
LCDcol = [0.2,0.2,0.2];
DPcol = [0.6,0.6,0.6];

%% plot fig7ai - CRT boxplot Rod

fig = figure('defaultAxesFontSize',12);
plotBoxplots(CRT,4,'PSDM_R (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,1.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 1];
print(fig, '..\supplementary_plots/figS2ai.pdf','-dpdf');

%% plot fig7aii - DP Rod boxplot

fig = figure('defaultAxesFontSize',12);
plotBoxplots(DP,4,'PSDM_R (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,1.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 1];
print(fig, '..\supplementary_plots/figS2aii.pdf','-dpdf');

%% plot fig7bi - CRT Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(CRT,Sim,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\supplementary_plots/figS2bi.pdf','-dpdf');

%% plot fig7bii - DP Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(DP,Sim,4,'PSDM_R (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\supplementary_plots/figS2bii.pdf','-dpdf');

%% plot fig7ci - CRT boxplot mel

fig = figure('defaultAxesFontSize',12);
plotBoxplots(CRT,5,'PSDM_I (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,1.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 1];
print(fig, '..\supplementary_plots/figS2ci.pdf','-dpdf');

%% plot fig7cii - DP mel boxplot

fig = figure('defaultAxesFontSize',12);
plotBoxplots(DP,5,'PSDM_I (%)')
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,1.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 1];
print(fig, '..\supplementary_plots/figS2cii.pdf','-dpdf');

%% plot fig7di - CRT Rod

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(CRT,Sim,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\supplementary_plots/figS2di.pdf','-dpdf');

%% plot fig7dii - DP Mel

fig = figure('defaultAxesFontSize',12);
h = plotDistortionsOverxy(DP,Sim,5,'PSDM_I (%)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\supplementary_plots/figS2dii.pdf','-dpdf');
%%
clear all;

%% functions

% plot spread
function h = plotDistortions(display,d,Sim,col,xlab,ylab)
    h=plot(Sim.ss(d,display.ssReproducible),display.ssDistorted(d,display.ssReproducible),'.','Color',col);
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
h=boxplot(display.distortionMetric(d,display.ssReproducible),'Orientation','horizontal','Colors','k','Symbol','k.');
set(h,{'linew'},{2});
grid on;
yticklabels({});
xlim([-100,100]);
xlabel(lab);
ylim([0.9,1.1]);
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
        xyDistortions(i,j)=mean(display.distortionMetric(d, display.ssReproducible & Sim.xyY(1,:)>x(i) & Sim.xyY(1,:)<x(i)+xStep & Sim.xyY(2,:)>y(j) & Sim.xyY(2,:)<y(j)+yStep));
    end
end
xyDistortions(isnan(xyDistortions)==1)=0;
h=imagesc(xyDistortions');
hold on;
plot(0.3128*100,0.3290*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('D10'));
c=colorbar;
c.Label.String = lab;
xlabel('CIE x');
ylabel('CIE y');
caxis([-60,60]);
xticks(1:10:100);
xticklabels(round(0:0.1:1,1));
yticks(1:10:100);
yticklabels(0:0.1:1);
xlim([10,70]);
ylim([10,70]);
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
        mbDistortions(i,j)=mean(display.distortionMetric(d,display.ssReproducible & display.mbDistorted(2,:)>mbx(i) & display.mbDistorted(2,:)<(mbx(i)+sStep) & display.mbDistorted(1,:)>mby(j) & display.mbDistorted(1,:)<(mby(j)+lStep)),'omitnan'); 
    end
end

mbDistortions(isnan(mbDistortions))=0;
h=imagesc(mbDistortions');
hold on;
plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('D10'));
c = colorbar();
c.Label.String = lab;
c.LineWidth = 2;
zlabel('CIE Y');
xlabel('L/(L+M)');
ylabel('S/(L+M)');
caxis([-60,60]);
xticks(1:10:100);
xticklabels(0:0.1:1);
yticks(1:10:100);
yticklabels(0:0.1:1);
ylim([0,15]);
xlim([50,100]);
axis square;
grid on;
end