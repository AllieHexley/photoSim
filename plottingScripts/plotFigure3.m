% plot Figure 3 of paper
% plots extension of MacLeod Boynton space for spectral and daylight locus
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% set up spectral locus plots

% Get the photoreceptor spectral sensitivities
% S, M, L, Rod, Mel
ss = GetCIES026;
wlsCIES026 = (390:1:780)';
T_cies026 = ss(:,11:end);
T_cies026(isnan(T_cies026)) = 0;

% set up MacLeod-Boynton chromaticity coordinates
lScale = 0.69283932;
mScale = 0.34967567;
sScale = 0.05547858;
% scale factors from CVRL MacLeod & Boynton (1979) 10-deg chromaticity
% coordinates based on the Stockman & Sharpe (2000) cone fundamentals: http://www.cvrl.org/

mb026(2,:) = T_cies026(2,:)*mScale;
mb026(3,:) = T_cies026(3,:)*lScale;
mb026(1,:) = T_cies026(1,:)*sScale;

iScale = 1./max(T_cies026(5,:)./(mb026(2,:)+mb026(3,:))); % so I/L+M peaks at 1
mb026(5,:) = T_cies026(5,:)*iScale;

% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
% scale for spds with 5nm spacing
mb026_5nm = mb026(:,1:5:end);
% remove Nans
mb026_5nm(isnan(mb026_5nm)) = 0;
mb026(isnan(mb026))=0;

% pick ccts to produce daylight locus for
cct = [4000, 5000, 5500, 6000, 6500, 7000, 7500, 8000, 9000, 10000, 11000, 12000, 13000];
DL = zeros(81, length(cct));
DLmbCoords = zeros(5, length(cct));
DLmb = zeros(3, length(cct));
for i=1:length(cct)
    DL = GenerateCIEDay(cct(i), B_cieday);
    A = trapz(380:5:780, DL);
    normDL = DL./A;
    DLmbCoords(:,i) = mb026_5nm*normDL(3:end);
    DLmb(1,i) = DLmbCoords(1,i)./(DLmbCoords(2,i)+DLmbCoords(3,i));
    DLmb(2,i) = DLmbCoords(3,i)./(DLmbCoords(2,i)+DLmbCoords(3,i));
    DLmb(3,i) = DLmbCoords(5,i)./(DLmbCoords(2,i)+DLmbCoords(3,i));
end

%% plot fig3a - Smb vs Lmb

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(SL.mb(2,:),SL.mb(1,:),'Color','k','LineWidth',2);

for i = [11,36,61,97,111,136,161,197,211,236,261,297,311,336,361] % mark on spectral locus crosses for each 25nm step
    plot(SL.mb(2,i),SL.mb(1,i),'x','Color',map(round(256*i/391),:),'MarkerSize',8,'LineWidth',2);
end

plot(DLmb(2,:),DLmb(1,:),'Color',[0,0,0.8],'LineWidth',2)
for i=1:2:13 % mark on crosses on daylight locus at every other CCT
    plot(DLmb(2,i),DLmb(1,i),'.','Color',map2(round(256*i/13),:),'LineWidth',2,'MarkerSize',12);
end
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
for i = [11,36,61,97,111,136,161,197,211,236,261,297,311,336,361] % mark on spectral locus crosses for each 25nm step
    plot(SL.mb(2,i),SL.mb(3,i),'x','Color',map(round(256*i/391),:),'MarkerSize',8,'LineWidth',2);
end
plot(DLmb(2,:),DLmb(3,:),'Color',[0,0,0.8],'LineWidth',2)
for i=1:2:13 % mark crosses on daylight locus
    plot(DLmb(2,i),DLmb(3,i),'.','Color',map2(round(256*i/13),:),'LineWidth',2,'MarkerSize',12);
end
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
for i = [11,36,61,97,111,136,161,197,211,236,261,297,311,336,361] % mark on spectral locus crosses for each 25nm step
    plot(SL.mb(3,i),SL.mb(1,i),'x','Color',map(round(256*i/391),:),'MarkerSize',8,'LineWidth',2);
end
plot(DLmb(3,:),DLmb(1,:),'Color',[0,0,0.8],'LineWidth',2)
for i=1:2:13 % plot crosses on daylight locus
    plot(DLmb(3,i),DLmb(1,i),'.','Color',map2(round(256*i/13),:),'LineWidth',2,'MarkerSize',12);
end
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