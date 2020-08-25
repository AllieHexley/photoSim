% script to plot RGB gamut for 3 primary display - RGB mythbusting
% created by ACH 01/07/2020

%% load data
load('photosimGamuts.mat');

%%  plot on xy diagram
figure('defaultAxesFontSize',18)
plotChromaticity();
hold on;
h(1)=plot(CRT.xyYMax(1,CRT.idx),CRT.xyYMax(2,CRT.idx),'k-','LineWidth',2);
h(2)=plot(LCD.xyYMax(1,LCD.idx),LCD.xyYMax(2,LCD.idx),'b-','LineWidth',2);
h(3) = plot(DP.xyYMax(1,DP.idx),DP.xyYMax(2,DP.idx),'-','Color',[0.75,0.75,0.75],'LineWidth',2);
legend(h,{'CRT','Display++','LCD'});
xlabel('x');
ylabel('y');
axis square

%% plot primaries out
figure('defaultAxesFontSize',18)
subplot(1,3,1)
g(1) = plot(390:780, CRT.rgb(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, CRT.rgb(:,2), 'g-','LineWidth',2);
plot(390:780, CRT.rgb(:,3), 'b-','LineWidth',2);
title('CRT');
xlabel('Wavelength (nm)');
ylabel('Radiance (uW/sr/m^2)');
xlim([390,780]);
ylim([0,0.011]);
axis square
subplot(1,3,2)
g(2) = plot(390:780, LCD.rgb(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, LCD.rgb(:,2), 'g-','LineWidth',2);
plot(390:780, LCD.rgb(:,3), 'b-','LineWidth',2);
title('LCD');
xlabel('Wavelength (nm)');
ylabel('Radiance (uW/sr/m^2)');
xlim([390,780]);
ylim([0,0.011]);
axis square
subplot(1,3,3)
g(3)=plot(390:780, DP.rgb(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, DP.rgb(:,2), 'g-','LineWidth',2);
plot(390:780, DP.rgb(:,3), 'b-','LineWidth',2);
title('Display++')
%legend(g,{'CRT','LCD','DP'});
xlabel('Wavelength (nm)');
ylabel('Radiance (uW/sr/m^2)');
xlim([390,780]);
ylim([0,0.011]);
axis square

%% plot five primaries
figure('defaultAxesFontSize',18)
subplot(1,2,1)
g(1) = plot(390:780, FP1.rgbcm(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, FP1.rgbcm(:,2), 'g-','LineWidth',2);
plot(390:780, FP1.rgbcm(:,3), 'b-','LineWidth',2);
plot(390:780, FP1.rgbcm(:,4), 'c-','LineWidth',2);
plot(390:780, FP1.rgbcm(:,5), 'm-','LineWidth',2);
title('Hypothetical Five Primary Display 1');
xlabel('Wavelength (nm)');
%ylabel('Radiance (uW/sr/m^2)');
xlim([390,780]);
yticklabels({});
%ylim([0,0.011]);
axis square

subplot(1,2,2)
g(1) = plot(390:780, FP2.rgbcm(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, FP2.rgbcm(:,2), 'g-','LineWidth',2);
plot(390:780, FP2.rgbcm(:,3), 'b-','LineWidth',2);
plot(390:780, FP2.rgbcm(:,4), 'c-','LineWidth',2);
plot(390:780, FP2.rgbcm(:,5), 'm-','LineWidth',2);
title('Hypothetical Five Primary Display 2');
xlabel('Wavelength (nm)');
%ylabel('Radiance (uW/sr/m^2)');
xlim([390,780]);
yticklabels({});
%ylim([0,0.011]);
axis square