% plot Figure 3 of paper
% plots extension of MacLeod Boynton space for spectral and daylight locus
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% plot fig3a - Smb vs Lmb

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(SL.mb(2,:),SL.mb(1,:),'Color','k','LineWidth',2);
xlabel('L/(L+M)');
ylabel('S/(L+M)');
%yticklabels({});
%xticks([400,500,600,700]);
%xticklabels({'400','500','600','700'});
%xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig3a.pdf','-dpdf');

%% plot fig3b - Imb vs Lmb

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(SL.mb(2,:),SL.mb(3,:),'Color','k','LineWidth',2);
xlabel('L/(L+M)');
ylabel('I/(L+M)');
%yticklabels({});
%xticks([400,500,600,700]);
%xticklabels({'400','500','600','700'});
%xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig3b.pdf','-dpdf');

%% plot fig3c - Smb vs Imb

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(SL.mb(3,:),SL.mb(1,:),'Color','k','LineWidth',2);
xlabel('I/(L+M)');
ylabel('S/(L+M)');
%yticklabels({});
%xticks([400,500,600,700]);
%xticklabels({'400','500','600','700'});
%xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig3c.pdf','-dpdf');

%%
clear all;