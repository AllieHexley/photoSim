% plot Figure 2 of paper
% plots correlations between photoreceptor signals in real world and
% correlation matrix
% created by ACH 01/07/2020

%% load data
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimMetrics_ReproduceLMS.mat');

%% plot fig2a - Correlation scatter plots

fig = figure('defaultAxesFontSize',12);
plotRealWorldCorrelations(Sim)
fig.PaperUnits = 'inches';
fig.PaperSize = [4.1,4.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 4 4];
print(fig, '..\plots\fig2a.png','-dpng');

%% plot fig2b -  Correlation matrix

fig = figure('defaultAxesFontSize',12);
plotRealWorldCorrelationHeatmaps(Sim)
fig.PaperUnits = 'inches';
fig.PaperSize = [4.1,4.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 4 4];
print(fig, '..\plots\fig2b.pdf','-dpdf');

%%

clear all

%% functions

function plotRealWorldCorrelations(Sim)

% set everything in lower diagonal to 0
i2K = [2,3,4,5,8,9,10,14,15,20];% indexes to keep
hold on;
c=1;
for i=1:5
    for j=1:5
        if ismember(c,i2K)==1
            subplot(5,5,c)
            plot(Sim.ss(i,:),Sim.ss(j,:),'k.');
            yticklabels({});
            xticklabels({});
            axis square
            box on;
            grid on;
        end
        c=c+1;
    end
end


end

function plotRealWorldCorrelationHeatmaps(Sim)

for ii=1:5
    % for real world
    [rhoS, pvalS] = corrcoef(Sim.ss(1,:),Sim.ss(ii,:));
    sCorr(ii) = rhoS(1,2);
        [rhoM, pvalM] = corrcoef(Sim.ss(2,:),Sim.ss(ii,:));
    mCorr(ii) = rhoM(1,2);
        [rhoL, pvalL] = corrcoef(Sim.ss(3,:),Sim.ss(ii,:));
    lCorr(ii) = rhoL(1,2);
        [rhoR, pvalR] = corrcoef(Sim.ss(4,:),Sim.ss(ii,:));
    rCorr(ii) = rhoR(1,2);
        [rhoI, pvalI] = corrcoef(Sim.ss(5,:),Sim.ss(ii,:));
    iCorr(ii) = rhoI(1,2);
end
% create matrix of correlations
correlationsTable = [sCorr;mCorr;lCorr;rCorr;iCorr];
% set everything in lower diagonal to 0
i2K = [6,11,12,16,17,18,21,22,23,24];% indexes to keep
for i=1:25
    if ismember(i,i2K)==0
        correlationsTable(i)=0;
    end
end
imagesc(correlationsTable);
hold on;
c=1;
for i=1:5
    for j=1:5
        if ismember(c,i2K)==1
            text(i-0.2,j,[num2str(round(correlationsTable(j,i),2))],'FontSize',10);
        end
        c=c+1;
    end
end
axis square
myColorMap = parula(256);
%myColorMap(1,:) = 1;
colormap(colorcet('gray'));%myColorMap);
hcb = colorbar;
ylabel(hcb,'\rho','FontSize',12,'Rotation',270);
hcb.Label.Position(1) = 3.5;
xticks([1:5]);
xticklabels(['S';'M';'L';'R';'I']);
yticks([1:5]);
yticklabels(['S';'M';'L';'R';'I']);

end