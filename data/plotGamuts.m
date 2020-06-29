% script to plot some colour gamut visualizations
% created by ACH on 16/06/2020

clear all;
close all;
clc;

%% load data

load('colorimetryNormPrimaries.mat');
load('spectralSensitivities.mat');

%% plot luminance correlation
figure('defaultAxesFontSize',18)

subplot(1,2,1)
plot(Sim.ss(3,:)+Sim.ss(2,:),Sim.ss(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('L+M (Photopic Lum)');
ylabel('R (Scotopic Lum)');
[rho, pval] = corrcoef(Sim.ss(3,:)+Sim.ss(2,:),Sim.ss(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(1,2,2)
plot(Sim.ss(3,:)+Sim.ss(2,:),Sim.ss(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('L+M (Photopic Lum)');
ylabel('I (Mel Irrad)');
[rho, pval] = corrcoef(Sim.ss(3,:)+Sim.ss(2,:),Sim.ss(5,:));
axis square
legend(num2str(rho(2,1)),'Location','best');

%% plot photoreceptor correlations

figure('Name','PRC','defaultAxesFontSize',18);
subplot(2,5,1)
plot(Sim.ss(1,:),Sim.ss(2,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('S');
ylabel('M');
[rho, pval] = corrcoef(Sim.ss(1,:),Sim.ss(2,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,2)
plot(Sim.ss(1,:),Sim.ss(3,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('S');
ylabel('L');
[rho, pval] = corrcoef(Sim.ss(1,:),Sim.ss(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,3)
plot(Sim.ss(1,:),Sim.ss(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('S');
ylabel('R');
[rho, pval] = corrcoef(Sim.ss(1,:),Sim.ss(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,4)
plot(Sim.ss(1,:),Sim.ss(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('S');
ylabel('I');
[rho, pval] = corrcoef(Sim.ss(1,:),Sim.ss(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,5)
plot(Sim.ss(2,:),Sim.ss(3,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('M');
ylabel('L');
[rho, pval] = corrcoef(Sim.ss(2,:),Sim.ss(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,6)
plot(Sim.ss(2,:),Sim.ss(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('M');
ylabel('R');
[rho, pval] = corrcoef(Sim.ss(2,:),Sim.ss(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,7)
plot(Sim.ss(2,:),Sim.ss(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('M');
ylabel('I');
[rho, pval] = corrcoef(Sim.ss(2,:),Sim.ss(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,8)
plot(Sim.ss(3,:),Sim.ss(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('L');
ylabel('R');
[rho, pval] = corrcoef(Sim.ss(3,:),Sim.ss(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,9)
plot(Sim.ss(3,:),Sim.ss(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('L');
ylabel('I');
[rho, pval] = corrcoef(Sim.ss(3,:),Sim.ss(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,10)
plot(Sim.ss(4,:),Sim.ss(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
xlabel('R');
ylabel('I');
[rho, pval] = corrcoef(Sim.ss(4,:),Sim.ss(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square


%% normalise photoreceptor responses to /L+M+S+R+I
norm = 0;
for i=1:5
    norm = norm+Sim.ss(i,:);
    Sim.ssNorm = Sim.ss./norm;
end

%% plot normalized photoreceptor correlations

figure('defaultAxesFontSize',18);
subplot(2,5,1)
plot(Sim.ssNorm(1,:),Sim.ssNorm(2,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('S/L+M+S+R+I');
ylabel('M/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(1,:),Sim.ssNorm(2,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,2)
plot(Sim.ssNorm(1,:),Sim.ssNorm(3,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('S/L+M+S+R+I');
ylabel('L/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(1,:),Sim.ssNorm(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,3)
plot(Sim.ssNorm(1,:),Sim.ssNorm(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('S/L+M+S+R+I');
ylabel('R/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(1,:),Sim.ssNorm(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,4)
plot(Sim.ssNorm(1,:),Sim.ssNorm(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('S/L+M+S+R+I');
ylabel('I/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(1,:),Sim.ssNorm(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,5)
plot(Sim.ssNorm(2,:),Sim.ssNorm(3,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('M/L+M+S+R+I');
ylabel('L/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(2,:),Sim.ssNorm(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,6)
plot(Sim.ssNorm(2,:),Sim.ssNorm(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('M/L+M+S+R+I');
ylabel('R/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(2,:),Sim.ssNorm(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,7)
plot(Sim.ssNorm(2,:),Sim.ssNorm(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('M/L+M+S+R+I');
ylabel('I/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(2,:),Sim.ssNorm(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,8)
plot(Sim.ssNorm(3,:),Sim.ssNorm(4,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('L/L+M+S+R+I');
ylabel('R/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(3,:),Sim.ssNorm(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,9)
plot(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('L/L+M+S+R+I');
ylabel('I/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(3,:),Sim.ssNorm(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,10)
plot(Sim.ssNorm(4,:),Sim.ssNorm(5,:),'kx','MarkerSize',0.02,'LineWidth',0.02);
xlabel('R/L+M+S+R+I');
ylabel('I/L+M+S+R+I');
[rho, pval] = corrcoef(Sim.ssNorm(4,:),Sim.ssNorm(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square

%% plot just 3 extended MB diagrams

figure('defaultAxesFontSize',18)
subplot(1,3,1)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
hold on;
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')
sgtitle('Unconstrained Simulated Spectra');

subplot(1,3,2)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
hold on;
xlabel('L/L+M+S+R+I')
ylabel('I/L+M+S+R+I')

subplot(1,3,3)
scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
hold on;
xlabel('I/L+M+S+R+I')
ylabel('S/L+M+S+R+I')

%% add in constraint that simulated spectra must be reproducible with display primaries
% i.e. go with constraint that we reproduce cones with priority, and see
% how these changes mel space

% calculate RGB to LMS conversion for the CRT display
% Calculate LMS cone values for each primary in the CRT (i.e. RGB to LMS
% matrix)
crtRGB2LMS = (T_cies026(1:3,:)*CRT.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
crtLMS2RGB = inv(crtRGB2LMS);

% calculate RGB to LMS conversion for the CRT display
% Calculate LMS cone values for each primary in the DP (i.e. RGB to LMS
% matrix)
DPRGB2LMS = (T_cies026(1:3,:)*DP.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
DPLMS2RGB = inv(DPRGB2LMS);

% calculate RGB to LMS conversion for the LCD display
% Calculate LMS cone values for each primary in the LCD (i.e. RGB to LMS
% matrix)
LCDRGB2LMS = (T_cies026(1:3,:)*LCD.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
LCDLMS2RGB = inv(LCDRGB2LMS);

%% for all simulated LMS values, calculate mel value and do quiver plots for how much mel value shifts
T_xyz = csvread('lin2012xyz10e_1_7sf.csv');
wls_xyz = T_xyz(:, 1);
T_xyz = 683*T_xyz(wls_xyz >= 390 & wls_xyz <= 780, 2:end)';
wls_xyz = wls_xyz(wls_xyz >= 390 & wls_xyz <= 780, 1);
wls_xyz = [];

for i = 1:length(Sim.ss)
    rgbSettingCRT = crtLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedCRT = (rgbSettingCRT'.*CRT.rgb);
    Sim.ssRealizedCRT(:,i) = T_cies026*(rgbRealizedCRT(:,1)+rgbRealizedCRT(:,2)+rgbRealizedCRT(:,3));
    Sim.xyRealizedCRT(:,i) = XYZToxyY(T_xyz*(rgbRealizedCRT(:,1)+rgbRealizedCRT(:,2)+rgbRealizedCRT(:,3)));
    
    rgbSettingDP = DPLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedDP = (rgbSettingDP'.*DP.rgb);
    Sim.ssRealizedDP(:,i) = T_cies026*(rgbRealizedDP(:,1)+rgbRealizedDP(:,2)+rgbRealizedDP(:,3));
    
    rgbSettingLCD = LCDLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedLCD = (rgbSettingLCD'.*LCD.rgb);
    Sim.ssRealizedLCD(:,i) = T_cies026*(rgbRealizedLCD(:,1)+rgbRealizedLCD(:,2)+rgbRealizedLCD(:,3));
end

%% plot on xy heatmap of how much did mel shift 

figure('defaultAxesFontSize',18)
%scatter(Sim.xyRealizedCRT(1,:),Sim.xyRealizedCRT(2,:),5,(Sim.ssRealizedCRT(5,:)-Sim.ss(5,:))./Sim.ss(5,:)),'filled');
scatter3(Sim.xyRealizedCRT(1,:),Sim.xyRealizedCRT(2,:),((Sim.ssRealizedCRT(5,:)-Sim.ss(5,:))./Sim.ss(5,:)).*100);
%plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
%c = colorbar;
%c.Label.String = 'Mel Actual - Mel Realized';
xlim([0,1]);
ylim([0,1]);
title('CRT');
xlabel('x');
ylabel('y');
zlabel('Mel Warp Contrast %');

%%

figure('defaultAxesFontSize',18)
scatter(Sim.xyRealizedCRT(1,:),Sim.xyRealizedCRT(2,:),5,100.*((Sim.ssRealizedCRT(5,:)-Sim.ss(5,:))./Sim.ss(5,:)),'filled');
%scatter3(Sim.xyRealizedCRT(1,:),Sim.xyRealizedCRT(2,:),(Sim.ssRealizedCRT(5,:)-Sim.ss(5,:))./Sim.ss(5,:));
%plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel Warp Contrast';
xlim([0,1]);
ylim([0,1]);
title('CRT');
xlabel('x');
ylabel('y');
%zlabel('Mel Warp Contrast');

%% plot warped spectra on correlations plot

figure(2);
subplot(2,5,1)
plot(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(2,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(1,idxCRT),CRT.ss(2,idxCRT),'g-','LineWidth',2);

xlabel('S');
ylabel('M');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(2,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,2)
plot(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(3,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(1,idxCRT),CRT.ss(3,idxCRT),'g-','LineWidth',2);

xlabel('S');
ylabel('L');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,3)
plot(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(4,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(1,idxCRT),CRT.ss(4,idxCRT),'g-','LineWidth',2);

xlabel('S');
ylabel('R');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,4)
plot(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(5,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(1,idxCRT),CRT.ss(5,idxCRT),'g-','LineWidth',2);

xlabel('S');
ylabel('I');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(1,:),Sim.ssRealizedCRT(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,5)
plot(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(3,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(2,idxCRT),CRT.ss(3,idxCRT),'g-','LineWidth',2);

xlabel('M');
ylabel('L');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(3,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,6)
plot(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(4,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(2,idxCRT),CRT.ss(4,idxCRT),'g-','LineWidth',2);

xlabel('M');
ylabel('R');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,7)
plot(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(5,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(2,idxCRT),CRT.ss(5,idxCRT),'g-','LineWidth',2);

xlabel('M');
ylabel('I');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(2,:),Sim.ssRealizedCRT(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,8)
plot(Sim.ssRealizedCRT(3,:),Sim.ssRealizedCRT(4,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(3,idxCRT),CRT.ss(4,idxCRT),'g-','LineWidth',2);

xlabel('L');
ylabel('R');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(3,:),Sim.ssRealizedCRT(4,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,9)
plot(Sim.ssRealizedCRT(3,:),Sim.ssRealizedCRT(5,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(3,idxCRT),CRT.ss(5,idxCRT),'g-','LineWidth',2);

xlabel('L');
ylabel('I');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(3,:),Sim.ssRealizedCRT(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square
subplot(2,5,10)
plot(Sim.ssRealizedCRT(4,:),Sim.ssRealizedCRT(5,:),'rx','MarkerSize',0.02,'LineWidth',0.02);
hold on;
%plot(CRT.ss(4,idxCRT),CRT.ss(5,idxCRT),'g-','LineWidth',2);

xlabel('R');
ylabel('I');
[rho, pval] = corrcoef(Sim.ssRealizedCRT(4,:),Sim.ssRealizedCRT(5,:));
legend(num2str(rho(2,1)),'Location','best');
axis square

%% plot input vs realized for each photoreceptor

figure('defaultAxesFontSize',18);
labels = {'S','M','L','R','I'};
mdl = cell(5,1);
for i =1:5
    subplot(3,5,i)
    plot(Sim.ss(i,:),Sim.ssRealizedCRT(i,:),'kx');
    xlabel('Desired');
    ylabel('Realized CRT');
    title(labels{i});
    axis square
    mdl{i} = fitlm(Sim.ss(i,:),Sim.ssRealizedCRT(i,:));
    text(0.06,0.01,['Rsq: ', num2str(round(mdl{i}.Rsquared.Ordinary,3))],'FontSize',18);
    xlim([0,0.12]);
    ylim([0,0.12]);
end

for i =1:5
    subplot(3,5,i+5)
    plot(Sim.ss(i,:),Sim.ssRealizedLCD(i,:),'bx');
    xlabel('Desired');
    ylabel('Realized LCD');
    title(labels{i});
    axis square
    mdl{i} = fitlm(Sim.ss(i,:),Sim.ssRealizedLCD(i,:));
    text(0.06,0.01,['Rsq: ', num2str(round(mdl{i}.Rsquared.Ordinary,3))],'FontSize',18);
    xlim([0,0.12]);
    ylim([0,0.12]);
end

for i =1:5
    subplot(3,5,10+i)
    plot(Sim.ss(i,:),Sim.ssRealizedDP(i,:),'rx');
    xlabel('Desired');
    ylabel('Realized DP');
    title(labels{i});
    axis square
    mdl{i} = fitlm(Sim.ss(i,:),Sim.ssRealizedDP(i,:));
    text(0.06,0.01,['Rsq: ', num2str(round(mdl{i}.Rsquared.Ordinary,3))],'FontSize',18);
    xlim([0,0.12]);
    ylim([0,0.12]);
end


%% calculate gamut in SML space of CRT

figure('defaultAxesFontSize',18)

hold on;
pairs = [1,2;1,3;1,4;1,5;,2,3;2,4;2,5;3,4;3,5;4,5];
pairNames = ['S','M';'S','L';'S','R';'S','I';'M','L';'M','R';'M','I';'L','R';'L','I';'R','I'];

for i = 1:10
    idxCRTcones{i}=convhull(CRT.ssFullRange(pairs(i,1),:),CRT.ssFullRange(pairs(i,2),:));
    idxLCDcones{i}=convhull(LCD.ssFullRange(pairs(i,1),:),LCD.ssFullRange(pairs(i,2),:));
    idxDPcones{i}=convhull(DP.ssFullRange(pairs(i,1),:),DP.ssFullRange(pairs(i,2),:));
end

totalSpec = 39699;
sgtitle('CRT');
for i=1:10
    subplot(2,5,i)
    p(1)=plot(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),'k.');
    hold on;
    p(2)=plot(Sim.ssRealizedCRT(pairs(i,1),:),Sim.ssRealizedCRT(pairs(i,2),:),'bx');
    %[rho{i}, pval{i}] = corrcoef(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:));
    %[rhoW{i}, pvalW{i}] = corrcoef(Sim.ssRealizedCRT(pairs(i,1),:),Sim.ssRealizedCRT(pairs(i,2),:));
    p(3)=plot(CRT.ssFullRange(pairs(i,1),idxCRTcones{i}),CRT.ssFullRange(pairs(i,2),idxCRTcones{i}),'g-','LineWidth',2);
    xlabel(pairNames(i,1));
    ylabel(pairNames(i,2));

    %legend(p,{'Sim','Realized'});
    
    [in{i},on{i}] = inpolygon(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),CRT.ssFullRange(pairs(i,1),idxCRTcones{i}),CRT.ssFullRange(pairs(i,2),idxCRTcones{i}));
    captured{i} = sum(in{i})+sum(on{i});
    percCaptured{i} = round((captured{i}./totalSpec).*100,1);
    
    [inWarp{i},onWarp{i}] = inpolygon(Sim.ssRealizedCRT(pairs(i,1),:),Sim.ssRealizedCRT(pairs(i,2),:),CRT.ssFullRange(pairs(i,1),idxCRTcones{i}),CRT.ssFullRange(pairs(i,2),idxCRTcones{i}));
    capturedWarp{i} = sum(inWarp{i})+sum(onWarp{i});
    percCapturedWarp{i} = round((capturedWarp{i}./totalSpec).*100,1);
    
    text(0.01,0.4,['% all sim: ', num2str(percCaptured{i})],'FontSize',18);
    text(0.01,0.45,['% warped: ', num2str(percCapturedWarp{i})],'FontSize',18);
    xlim([0,0.5]);
    ylim([0,0.5]);
end

%% and for LCD
figure('defaultAxesFontSize',18)
sgtitle('LCD');
for i=1:10
    subplot(2,5,i)
    p(1)=plot(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),'k.');
    hold on;
    p(2)=plot(Sim.ssRealizedLCD(pairs(i,1),:),Sim.ssRealizedLCD(pairs(i,2),:),'bx');
    p(3)=plot(LCD.ssFullRange(pairs(i,1),idxLCDcones{i}),LCD.ssFullRange(pairs(i,2),idxLCDcones{i}),'g-','LineWidth',2);
    xlabel(pairNames(i,1));
    ylabel(pairNames(i,2));

    %legend(p,{'Sim','Realized','LCD gamut'});
    
    [in{i},on{i}] = inpolygon(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),LCD.ssFullRange(pairs(i,1),idxLCDcones{i}),LCD.ssFullRange(pairs(i,2),idxLCDcones{i}));
    captured{i} = sum(in{i})+sum(on{i});
    percCaptured{i} = round((captured{i}./totalSpec).*100,1);
    
    [inWarp{i},onWarp{i}] = inpolygon(Sim.ssRealizedLCD(pairs(i,1),:),Sim.ssRealizedLCD(pairs(i,2),:),LCD.ssFullRange(pairs(i,1),idxLCDcones{i}),LCD.ssFullRange(pairs(i,2),idxLCDcones{i}));
    capturedWarp{i} = sum(inWarp{i})+sum(onWarp{i});
    percCapturedWarp{i} = round((capturedWarp{i}./totalSpec).*100,1);
    
    text(0.01,0.4,['% all sim: ', num2str(percCaptured{i})],'FontSize',18);
    text(0.01,0.45,['% warped: ', num2str(percCapturedWarp{i})],'FontSize',18);
    xlim([0,0.5]);
    ylim([0,0.5]);
end

%% and DP
figure('defaultAxesFontSize',18)
sgtitle('Display ++');
for i=1:10
    subplot(2,5,i)
    p(1)=plot(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),'k.');
    hold on;
    p(2)=plot(Sim.ssRealizedDP(pairs(i,1),:),Sim.ssRealizedDP(pairs(i,2),:),'bx');
    p(3)=plot(DP.ssFullRange(pairs(i,1),idxDPcones{i}),DP.ssFullRange(pairs(i,2),idxDPcones{i}),'g-','LineWidth',2);
    xlabel(pairNames(i,1));
    ylabel(pairNames(i,2));

    %legend(p,{'Sim','Realized','DP gamut'});
    
    [in{i},on{i}] = inpolygon(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),DP.ssFullRange(pairs(i,1),idxDPcones{i}),DP.ssFullRange(pairs(i,2),idxDPcones{i}));
    captured{i} = sum(in{i})+sum(on{i});
    percCaptured{i} = round((captured{i}./totalSpec).*100,1);
    
    [inWarp{i},onWarp{i}] = inpolygon(Sim.ssRealizedDP(pairs(i,1),:),Sim.ssRealizedDP(pairs(i,2),:),DP.ssFullRange(pairs(i,1),idxDPcones{i}),DP.ssFullRange(pairs(i,2),idxDPcones{i}));
    capturedWarp{i} = sum(inWarp{i})+sum(onWarp{i});
    percCapturedWarp{i} = round((capturedWarp{i}./totalSpec).*100,1);
    
    text(0.01,0.4,['% all sim: ', num2str(percCaptured{i})],'FontSize',18);
    text(0.01,0.45,['% warped: ', num2str(percCapturedWarp{i})],'FontSize',18);
    xlim([0,0.5]);
    ylim([0,0.5]);
end


%% normalise to find /L+M+S+R+I

normCRT = 0;
normDP = 0;
normLCD = 0;
normSim = 0;
normCRTFR = 0;
for i=1:5
    normCRT = normCRT+Sim.ssRealizedCRT(i,:);
    Sim.ssRealizedNormCRT = Sim.ssRealizedCRT./normCRT;
    normSim = normSim+Sim.ss(i,:);
    Sim.ssRealizedCRTNormSim = Sim.ssRealizedCRT./normSim;
    normCRTFR = normCRTFR+CRT.ssFullRange(i,:);
    CRT.ssFullRangeNorm = CRT.ssFullRange./normCRTFR;
    normDP = normDP+Sim.ssRealizedDP(i,:);
    Sim.ssRealizedNormDP = Sim.ssRealizedDP./normDP;
    normLCD = normLCD+Sim.ssRealizedLCD(i,:);
    Sim.ssRealizedNormLCD = Sim.ssRealizedLCD./normLCD;
end

%% plot both types of normalization - add correlations to this plot! - not sure what meaningful gamut is here
clear rho pval
CRT.ssFullRangeNorm(isnan(CRT.ssFullRangeNorm))=0;
for i = 1:10
    idxCRTconesnorm{i}=convhull(CRT.ssFullRangeNorm(pairs(i,1),:),CRT.ssFullRangeNorm(pairs(i,2),:));
end
figure('defaultAxesFontSize',18)
sgtitle('CRT');
for i=1:10
    subplot(2,5,i)
    p(1)=plot(Sim.ssNorm(pairs(i,1),:),Sim.ssNorm(pairs(i,2),:),'k.');
    hold on;
    p(2)=plot(Sim.ssRealizedCRTNormSim(pairs(i,1),:),Sim.ssRealizedCRTNormSim(pairs(i,2),:),'bx');
    [rho{i}, pval{i}] = corrcoef(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:));
    [rhoW{i}, pvalW{i}] = corrcoef(Sim.ssRealizedCRT(pairs(i,1),:),Sim.ssRealizedCRT(pairs(i,2),:));
    %p(3)=plot(CRT.ssFullRangeNorm(pairs(i,1),idxCRTcones{i}),CRT.ssFullRangeNorm(pairs(i,2),idxCRTcones{i}),'g-','LineWidth',2);
    xlabel([pairNames(i,1),'./S+M+L+R+I']);
    ylabel([pairNames(i,2),'./S+M+L+R+I']);

    legend(p,{'Sim','Realized'});
    
%     [in{i},on{i}] = inpolygon(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:),CRT.ssFullRange(pairs(i,1),idxCRTcones{i}),CRT.ssFullRange(pairs(i,2),idxCRTcones{i}));
%     captured{i} = sum(in{i})+sum(on{i});
%     percCaptured{i} = round((captured{i}./totalSpec).*100,1);
%     
%     [inWarp{i},onWarp{i}] = inpolygon(Sim.ssRealizedCRT(pairs(i,1),:),Sim.ssRealizedCRT(pairs(i,2),:),CRT.ssFullRange(pairs(i,1),idxCRTcones{i}),CRT.ssFullRange(pairs(i,2),idxCRTcones{i}));
%     capturedWarp{i} = sum(inWarp{i})+sum(onWarp{i});
%     percCapturedWarp{i} = round((capturedWarp{i}./totalSpec).*100,1);
    
    text(0.3,0.6,['Sim: ', num2str(round(rho{1,i}(2),2))],'FontSize',18);
    text(0.3,0.45,['Warped: ', num2str(round(rhoW{1,i}(2),2))],'FontSize',18);
    xlim([0,0.7]);
    ylim([0,0.7]);

end

% %%
% figure('defaultAxesFontSize',18)
% subplot(1,3,1)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),'r.');
% idxCRT = convhull(CRT.xyY(1,:), CRT.xyY(2,:));
% %k(1)=plot(CRT.ss(3,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(1,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-','LineWidth',2);
% xlabel('L/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% sgtitle('Constrained Simulated Spectra');
% 
% subplot(1,3,2)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(5,:),'r.');
% %k(1)=plot(CRT.ss(3,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(5,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-','LineWidth',2);
% xlabel('L/L+M+S+R+I')
% ylabel('I/L+M+S+R+I')
% title('CRT')
% 
% subplot(1,3,3)
% scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormCRT(5,:),Sim.ssRealizedNormCRT(1,:),'r.');
% %k(1)=plot(CRT.ss(5,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(1,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-','LineWidth',2);
% xlabel('I/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% 
% %% calculate RMS between desired and resulting mel level and plot as histogram of RMS levels
% 
% % some error happening here...
% a=[Sim.ssRealizedCRT(1,:);Sim.ss(1,:)];
% arms = rms(a,1);
% b=[Sim.ssRealizedCRT(2,:);Sim.ss(2,:)];
% brms = rms(b,1);
% c=[Sim.ssRealizedCRT(3,:);Sim.ss(3,:)];
% crms = rms(c,1);
% d=[Sim.ssRealizedCRT(4,:);Sim.ss(4,:)];
% drms = rms(d,1);
% e=[Sim.ssRealizedCRT(5,:);Sim.ss(5,:)];
% erms = rms(e,1);
% figure()
% subplot(3,2,1)
% histogram(arms);
% subplot(3,2,2)
% histogram(brms);
% subplot(3,2,3)
% histogram(crms);
% subplot(3,2,4)
% histogram(drms);
% subplot(3,2,5)
% histogram(erms);
% 
% 
% % subplot(3,3,4)
% % scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormDP(3,:),Sim.ssRealizedNormDP(1,:),'r.');
% % xlabel('L/L+M+S+R+I')
% % ylabel('S/L+M+S+R+I')
% % sgtitle('Constrained Simulated Spectra');
% % 
% % subplot(3,3,5)
% % scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormDP(3,:),Sim.ssRealizedNormDP(5,:),'r.');
% % xlabel('L/L+M+S+R+I')
% % ylabel('I/L+M+S+R+I')
% % title('DP');
% % 
% % subplot(3,3,6)
% % scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormDP(5,:),Sim.ssRealizedNormDP(1,:),'r.');
% % xlabel('I/L+M+S+R+I')
% % ylabel('S/L+M+S+R+I')
% % 
% % subplot(3,3,7)
% % scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormLCD(3,:),Sim.ssRealizedNormLCD(1,:),'r.');
% % xlabel('L/L+M+S+R+I')
% % ylabel('S/L+M+S+R+I')
% % sgtitle('Constrained Simulated Spectra');
% % 
% % subplot(3,3,8)
% % scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormLCD(3,:),Sim.ssRealizedNormLCD(5,:),'r.');
% % xlabel('L/L+M+S+R+I')
% % ylabel('I/L+M+S+R+I')
% % title('LCD');
% % 
% % subplot(3,3,9)
% % scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% % hold on;
% % scatter(Sim.ssRealizedNormLCD(5,:),Sim.ssRealizedNormLCD(1,:),'r.');
% % xlabel('I/L+M+S+R+I')
% % ylabel('S/L+M+S+R+I')
% 
% %% plot as quiver plots
% 
% figure('defaultAxesFontSize',18)
% subplot(1,3,1)
% %scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% hold on;
% b=100;
% quiver(Sim.ssNorm(3,1:b:39600),Sim.ssNorm(1,1:b:39600),Sim.ssRealizedNormCRT(3,1:b:39600)-Sim.ssNorm(3,1:b:39600),Sim.ssRealizedNormCRT(1,1:b:39600)-Sim.ssNorm(1,1:b:39600),0,'MaxHeadSize',0.2);
% %scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),'r.','LineWidth',0.5);
% xlabel('L/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% sgtitle('Constrained Simulated Spectra');
% 
% subplot(1,3,2)
% b=800;
% quiver(Sim.ssNorm(3,1:b:39600),Sim.ssNorm(5,1:b:39600),Sim.ssRealizedNormCRT(3,1:b:39600)-Sim.ssNorm(3,1:b:39600),Sim.ssRealizedNormCRT(5,1:b:39600)-Sim.ssNorm(5,1:b:39600),0,'MaxHeadSize',0.2);
% xlabel('L/L+M+S+R+I')
% ylabel('I/L+M+S+R+I')
% title('CRT')
% 
% subplot(1,3,3)
% quiver(Sim.ssNorm(5,1:b:39600),Sim.ssNorm(1,1:b:39600),Sim.ssRealizedNormCRT(5,1:b:39600)-Sim.ssNorm(5,1:b:39600),Sim.ssRealizedNormCRT(1,1:b:39600)-Sim.ssNorm(1,1:b:39600),0,'MaxHeadSize',0.2);
% xlabel('I/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% 
% %% try to calculate RMS of both original and warped spaces and RMS of difference between points
% 
% % RMS of the simulated data points
% %rmsSimLS = rms(Sim.ssNorm(
% 
% %% plot as 3D surface plot
% 
% normCRTCombo = 0;
% for i=1:5
%     normCRTCombo = normCRTCombo+CRT.ssFullRange(i,:);
%     CRT.ssFullRangeNorm = CRT.ssFullRange./normCRTCombo;
% end
% CRT.ssFullRangeNorm(isnan(CRT.ssFullRangeNorm)==1)=0;
% 
% figure('defaultAxesFontSize',18)
% subplot(1,2,1)
% title('Original');
% scatter3(Sim.ssNorm(3,:),Sim.ssNorm(1,:),Sim.ssNorm(5,:),'k.');
% xlabel('L/S+M+L+R+I');
% ylabel('S/S+M+L+R+I');
% zlabel('I/S+M+L+R+I');
% hold on;
% subplot(1,2,2)
% title('Warped to Display Gamut');
% scatter3(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),Sim.ssRealizedNormCRT(5,:),'r.');
% xlabel('L/S+M+L+R+I');
% ylabel('S/S+M+L+R+I');
% zlabel('I/S+M+L+R+I');
% hold on;
% 
% % attempt to produce mesh
% [xsim,ysim] = meshgrid(0:0.01:0.8,0:0.01:0.4);
% zsim = griddata(Sim.ssNorm(3,:),Sim.ssNorm(1,:),Sim.ssNorm(5,:),xsim,ysim);
% subplot(1,2,1)
% M = mesh(xsim,ysim,zsim);
% 
% % attempt to produce mesh for Realized Data
% [xcrt,ycrt] = meshgrid(0:0.01:0.8,0:0.01:0.4);
% zcrt = griddata(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),Sim.ssRealizedNormCRT(5,:),xcrt,ycrt);
% subplot(1,2,2)
% mesh(xcrt,ycrt,zcrt);
% 
% %%
% 
% % attempt to produce mesh for gamut of CRT
% [xcrtCombo,ycrtCombo] = meshgrid(0:0.01:0.8,0:0.01:0.4);
% zcrtCombo = griddata(CRT.ssFullRangeNorm(3,:),CRT.ssFullRangeNorm(1,:),CRT.ssFullRangeNorm(5,:),xcrtCombo,ycrtCombo);
% %mesh(xcrtCombo,ycrtCombo,zcrtCombo);
% 
% shp = alphaShape(Sim.ssNorm(3,:)',Sim.ssNorm(1,:)',Sim.ssNorm(5,:)');
% volume(shp)
% shpCRT = alphaShape(Sim.ssRealizedNormCRT(3,:)',Sim.ssRealizedNormCRT(1,:)',Sim.ssRealizedNormCRT(5,:)');
% volume(shpCRT)
% % this doesn't seem to work the way I want it to...
% % perhaps try and look at Rafal's code for producing 3D gamut for MPHDR
% % display here
% shpCRTCombo = alphaShape(CRT.ssFullRangeNorm(3,:)',CRT.ssFullRangeNorm(1,:)',CRT.ssFullRangeNorm(5,:)');
% volume(shpCRTCombo)
% 
% % possibly plot on 2D projections in all possible combinations with spread
% % along all photoreceptor axes 
% 
% %% try with SLI matrix instead of LMS matrix
% 
% % calculate RGB to LMS conversion for the CRT display
% % Calculate LMS cone values for each primary in the CRT (i.e. RGB to LMS
% % matrix)
% crtRGB2SLI = (T_cies026([1,3,5],:)*CRT.rgb);
% 
% % find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
% crtLMS2RGBSLI = inv(crtRGB2SLI);
% 
% 
% %% do the same thing but with and SLI matrix instead of an LMS matrix
% 
% for i = 1:length(Sim.ss)
%     rgbSettingCRT = crtLMS2RGBSLI*Sim.ss([1,3,5],i);
%     rgbRealizedCRT = (rgbSettingCRT'.*CRT.rgb);
%     Sim.ssRealizedSLICRT(:,i) = T_cies026*(rgbRealizedCRT(:,1)+rgbRealizedCRT(:,2)+rgbRealizedCRT(:,3));
% end
% 
% normCRT = 0;
% for i=1:5
%     normCRT = normCRT+Sim.ssRealizedSLICRT(i,:);
%     Sim.ssRealizedNormSLICRT = Sim.ssRealizedSLICRT./normCRT;
% end
% 
% 
% figure('defaultAxesFontSize',18)
% subplot(1,3,1)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormSLICRT(3,:),Sim.ssRealizedNormSLICRT(1,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% sgtitle('Constrained Simulated Spectra');
% 
% subplot(1,3,2)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormSLICRT(3,:),Sim.ssRealizedNormSLICRT(5,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('I/L+M+S+R+I')
% title('CRT')
% 
% subplot(1,3,3)
% scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormSLICRT(5,:),Sim.ssRealizedNormSLICRT(1,:),'r.');
% xlabel('I/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')