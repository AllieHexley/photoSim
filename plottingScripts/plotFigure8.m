% plot Figure 9 of paper
% plots MB space and distorted MB
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% plot full real-world database spectra in MacLeod-Boynton diagram and two new projections

plotMBSpace(Sim, 'a')

%% plot distorted MacLeod-Boynton colour space for all five displays
plotDistortedColourMaps(CRT,Sim,'b_crt_')
plotDistortedColourMaps(LCD,Sim,'b_lcd_')
plotDistortedColourMaps(DP,Sim,'b_dp_')
plotDistortedColourMaps(nb5p,Sim,'b_nb5p_')
plotDistortedColourMaps(bb5p,Sim,'b_bb5p_')
%%
clear all;

%% functions

function plotMBSpace(Sim, fignum)

lStep = 0.003;
sStep = 0.001;
mbx = 0.6:0.003:0.9;
mby = 0:0.001:0.1;
for i=1:100
    for j=1:100
        mbFreq(i,j)=sum((Sim.mb(2,:)>mbx(i) & Sim.mb(2,:)<(mbx(i)+sStep) & Sim.mb(1,:)>mby(j) & Sim.mb(1,:)<(mby(j)+lStep)),'omitnan'); 
    end
end

mbFreq(isnan(mbFreq))=0;
%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(mbFreq');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
xlabel('L/(L+M)');
ylabel('S/(L+M)');
caxis([0,200]);
xticks([1,33,66,99]);
xticklabels([0.6,0.7,0.8,0.9]);
yticks([1,50,100]);
yticklabels([0,0.05,0.1]);
%ylim([0,30]);
%xlim([50,200]);
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
set(gca,'Color','b')
print(fig, ['..\plots\fig8' fignum 'i.pdf'],'-dpdf');

%% second projection
lStep = 0.003;
sStep = 0.0025;
mbx = 0.6:0.003:0.9;
mby = 0:0.0025:0.25;
for i=1:100
    for j=1:100
        mbFreq2(i,j)=sum((Sim.mb(2,:)>mbx(i) & Sim.mb(2,:)<(mbx(i)+sStep) & Sim.mb(3,:)>mby(j) & Sim.mb(3,:)<(mby(j)+lStep)),'omitnan'); 
    end
end

mbFreq2(isnan(mbFreq2))=0;

%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(mbFreq2');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
xlabel('L/(L+M)');
ylabel('I/(L+M)');
caxis([0,200]);
xticks([1,33,66,99]);
xticklabels([0.6,0.7,0.8,0.9]);
yticks([1,40,80]);
yticklabels([0,0.1,0.2]);
%ylim([0,30]);
%xlim([50,200]);
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, ['..\plots\fig8' fignum 'ii.pdf'],'-dpdf');

%% third projection
sStep = 0.001;
lStep = 0.0025;
mby = 0:0.001:0.1;
mbx = 0:0.0025:0.25;
for i=1:100
    for j=1:100
        mbFreq3(i,j)=sum((Sim.mb(3,:)>mbx(i) & Sim.mb(3,:)<(mbx(i)+sStep) & Sim.mb(1,:)>mby(j) & Sim.mb(1,:)<(mby(j)+lStep)),'omitnan'); 
    end
end

mbFreq3(isnan(mbFreq3))=0;

%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(mbFreq3');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
xlabel('I/(L+M)');
ylabel('S/(L+M)');
caxis([0,200]);
yticks([1,50,100]);
yticklabels([0,0.05,0.1]);
xticks([1,40,80]);
xticklabels([0,0.1,0.2]);
%ylim([0,30]);
%xlim([50,200]);
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, ['..\plots\fig8' fignum 'iii.pdf'],'-dpdf');

end

function plotDistortedColourMaps(disp,Sim, fignum)

ifReproducible = disp.ssReproducible;
% check if within 1% of error for each signal
withinTolerance = (disp.ssDistorted+(disp.ssDistorted*0.01)) >= Sim.ss & (disp.ssDistorted-(disp.ssDistorted*0.01)) <= Sim.ss; % to 1% tolerance
ifWithinTolerance = (sum(withinTolerance(:,:))==5);
% check if within tolerance and reproducible
ifMatch = ifWithinTolerance+ifReproducible;

%%
lStep = 0.003;
sStep = 0.001;
mbx = 0.6:0.003:0.9;
mby = 0:0.001:0.1;
for i=1:100
    for j=1:100
        disp1(i,j)=sum((disp.mbDistorted(2,ifMatch==2)>mbx(i) & disp.mbDistorted(2,ifMatch==2)<(mbx(i)+sStep) & disp.mbDistorted(1,ifMatch==2)>mby(j) & disp.mbDistorted(1,ifMatch==2)<(mby(j)+lStep)),'omitnan'); 
    end
end

disp1(isnan(disp1))=0;
%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(disp1');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
%c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
%xlabel('L/(L+M)');
%ylabel('S/(L+M)');
caxis([0,200]);
xticks([1,33,66,99]);
xticklabels([0.6,0.7,0.8,0.9]);
yticks([1,50,100]);
yticklabels([0,0.05,0.1]);
% ylim([0,30]);
% xlim([50,200]);
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, ['..\plots\fig8' fignum 'i.pdf'],'-dpdf');

%% second projection
lStep = 0.003;
sStep = 0.0025;
mbx = 0.6:0.003:0.9;
mby = 0:0.0025:0.25;
for i=1:100
    for j=1:100
        disp2(i,j)=sum((disp.mbDistorted(2,ifMatch==2)>mbx(i) & disp.mbDistorted(2,ifMatch==2)<(mbx(i)+sStep) & disp.mbDistorted(3,ifMatch==2)>mby(j) & disp.mbDistorted(3,ifMatch==2)<(mby(j)+lStep)),'omitnan'); 
    end
end

disp2(isnan(disp2))=0;

%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(disp2');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
%c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
%xlabel('L/(L+M)');
%ylabel('I/(L+M)');
caxis([0,200]);
xticks([1,33,66,99]);
xticklabels([0.6,0.7,0.8,0.9]);
yticks([1,40,80]);
yticklabels([0,0.1,0.2]);
% ylim([0,30]);
% xlim([50,200]);
%xticklabels({});
%yticklabels({});
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, ['..\plots\fig8' fignum 'ii.pdf'],'-dpdf');

%% third projection
sStep = 0.001;
lStep = 0.0025;
mby = 0:0.001:0.1;
mbx = 0:0.0025:0.25;
for i=1:100
    for j=1:100
        disp3(i,j)=sum((disp.mbDistorted(3,ifMatch==2)>mbx(i) & disp.mbDistorted(3,ifMatch==2)<(mbx(i)+sStep) & disp.mbDistorted(1,ifMatch==2)>mby(j) & disp.mbDistorted(1,ifMatch==2)<(mby(j)+lStep)),'omitnan'); 
    end
end

disp3(isnan(disp3))=0;

%%
fig = figure('defaultAxesFontSize',12);
h=imagesc(disp3');
hold on;
% plot(0.6922*100,0.0257*100,'ko','LineWidth',2,'MarkerSize',2);
set(gca,'YDir','normal');
colormap(colorcet('L17'));
%c = colorbar();
c.Label.String = 'Frequency';
c.LineWidth = 2;
zlabel('CIE Y');
%xlabel('I/(L+M)');
%ylabel('S/(L+M)');
caxis([0,200]);
yticks([1,50,100]);
yticklabels([0,0.05,0.1]);
xticks([1,40,80]);
xticklabels([0,0.1,0.2]);
% ylim([0,30]);
% xlim([50,200]);
axis square;
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, ['..\plots\fig8' fignum 'iii.pdf'],'-dpdf');

end



