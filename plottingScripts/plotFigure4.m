% plot Figure 4 of paper
% plots primaries of displays used
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% plot Fig3a - CRT primaries

fig = figure('defaultAxesFontSize',12);
plotDisplayPrimaries(CRT)
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig4a.pdf','-dpdf');

%% plot Fig3b -  Dell LCD primaries

fig = figure('defaultAxesFontSize',12);
plotDisplayPrimaries(LCD)
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig4b.pdf','-dpdf');

%% plot Fig3c - Display++ primaries

fig = figure('defaultAxesFontSize',12);
plotDisplayPrimaries(DP)
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig4c.pdf','-dpdf');

%% plot Fig3d - Broadband five primary primaries

fig = figure('defaultAxesFontSize',12);
plotDisplayPrimariesMultiPrimary(nb5p)
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig4d.pdf','-dpdf');

%% plot Fig3e - Narrowband five primary primaries

fig = figure('defaultAxesFontSize',12);
plotDisplayPrimariesMultiPrimary(bb5p)
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig4e.pdf','-dpdf');

%%
clear all;

%% functions

function plotDisplayPrimaries(display)

hold on;
h(1)=plot(390:780,display.spd(:,1),'Color',[0.5,0,0],'LineWidth',2);
h(2)=plot(390:780,display.spd(:,2),'Color',[0,0.5,0],'LineWidth',2);
h(3)=plot(390:780,display.spd(:,3),'Color',[0,0,0.5],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;

end

function plotDisplayPrimariesMultiPrimary(display)

hold on;
h(1)=plot(390:780,display.spd(:,5),'Color',[0.5,0,0],'LineWidth',2);
h(2)=plot(390:780,display.spd(:,3),'Color',[0,0.5,0],'LineWidth',2);
h(3)=plot(390:780,display.spd(:,2),'Color',[0,0,0.5],'LineWidth',2);
h(3)=plot(390:780,display.spd(:,1),'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h(3)=plot(390:780,display.spd(:,4),'Color',[0.9290,0.6940,0.1250],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;

end