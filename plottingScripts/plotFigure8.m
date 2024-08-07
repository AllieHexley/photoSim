% plot Supplementary Figure 2 of paper
% plots MB space and distorted MB
% created by ACH 13/07/2021

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('../photosimMetrics_ReproduceLMS.mat');

%% plot distorted MacLeod-Boynton colour space for the five primary displays

%plotDistortedColourMaps(MPHDR,Sim,'MPHDR')
%plotDistortedColourMaps(Man,Sim,'Man')
plotDistortedColourMaps(NZ,Sim,'NZ')

%%
clear all;

%% functions

function plotDistortedColourMaps(disp,Sim, fignum)

ifReproducible = disp.ssReproducible;
% check if within 1% of error for each signal
withinTolerance = (disp.ssDistorted+(disp.ssDistorted*0.01)) >= Sim.ss & (disp.ssDistorted-(disp.ssDistorted*0.01)) <= Sim.ss; % to 1% tolerance
ifWithinTolerance = (sum(withinTolerance(:,:))==5);
% check if within tolerance and reproducible
if size(ifWithinTolerance,1)<size(ifWithinTolerance,2)
    ifWithinTolerance = ifWithinTolerance';
end
if size(ifReproducible,1)<size(ifReproducible,2)
    ifReproducible = ifReproducible';
end
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
print(fig, ['..\plots/SF3' fignum 'i.pdf'],'-dpdf');

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
print(fig, ['..\plots/SF3' fignum 'ii.pdf'],'-dpdf');

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
print(fig, ['..\plots/SF3' fignum 'iii.pdf'],'-dpdf');

end



