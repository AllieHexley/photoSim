% plot Figure 5 of paper
% plots MB space and distorted MB
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
nb5pcol = [0, 0.8, 0];
bb5pcol = [0, 0, 0.8];

%% plot fig8ai - MB space real world

fig = figure('defaultAxesFontSize',12);
plot(Sim.mb(2,:),Sim.mb(1,:),'.','Color',[0.3,0,0.6]);
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8ai.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plot(Sim.mb(2,:),Sim.mb(3,:),'.','Color',[0.3,0,0.6]);
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8aii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plot(Sim.mb(3,:),Sim.mb(1,:),'.','Color',[0.3,0,0.6]);
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8aiii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(CRT,2,1,CRTcol,Sim)
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8bi.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(CRT,2,3,CRTcol,Sim)
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8bii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(CRT,3,1,CRTcol,Sim)
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8biii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(LCD,2,1,LCDcol,Sim)
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8ci.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(LCD,2,3,LCDcol,Sim)
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8cii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(LCD,3,1,LCDcol,Sim)
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8ciii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(DP,2,1,DPcol,Sim)
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8di.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(DP,2,3,DPcol,Sim)
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8dii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(DP,3,1,DPcol,Sim)
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8diii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(nb5p,2,1,nb5pcol,Sim)
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8ei.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(nb5p,2,3,nb5pcol,Sim)
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8eii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(nb5p,3,1,nb5pcol,Sim)
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8eiii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(bb5p,2,1,bb5pcol,Sim)
xlabel('L/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8fi.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(bb5p,2,3,bb5pcol,Sim)
xlabel('L/(L+M)');
ylabel('I/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8fii.pdf','-dpdf');

%% plot fig8aii - MB space real world

fig = figure('defaultAxesFontSize',12);
plotMBSpace(bb5p,3,1,bb5pcol,Sim)
xlabel('I/(L+M)');
ylabel('S/(L+M)');
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig8fiii.pdf','-dpdf');

%%
clear all;

%% functions

function plotMBSpace(display,x,y,col,Sim)

ifReproducible = display.ssReproducible;
% check if within 1% of error for each signal
withinTolerance = (display.ssDistorted+(display.ssDistorted*0.01)) >= Sim.ss & (display.ssDistorted-(display.ssDistorted*0.01)) <= Sim.ss; % to 1% tolerance
ifWithinTolerance = (sum(withinTolerance(:,:))==5);
% check if within tolerance and reproducible
ifMatch = ifWithinTolerance+ifReproducible;
plot(display.mbDistorted(x,ifMatch==2),display.mbDistorted(y,ifMatch==2),'.','Color',col);

end