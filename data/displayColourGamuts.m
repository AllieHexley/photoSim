% script to convert RGB primaries to xy and xyZ space and to
% MacLeod-Boynton space to see impact of melanopsin on gamut shape
% created by ACH on 13/05/2020

clear all;
clc;

%% Set up colorimetry
% Get the CIE 2015 10degree XYZ functions
T_xyz = csvread('lin2012xyz10e_1_7sf.csv');
wls_xyz = T_xyz(:, 1);
T_xyz = 683*T_xyz(wls_xyz >= 390 & wls_xyz <= 780, 2:end)';
wls_xyz = wls_xyz(wls_xyz >= 390 & wls_xyz <= 780, 1);
wls_xyz = [];

% Get the photoreceptor spectral sensitivities
% S, M, L, Rod, Mel
ss = GetCIES026;
wlsCIES026 = (390:1:780)';
T_cies026 = ss(:,11:end);
T_cies026(isnan(T_cies026)) = 0;

%% load in the primaries for a particular display (let's start with the CRT)
load('RGBPhospher.mat');
wlsCRT = RGBPhospher(:,1);
rgbCRT = RGBPhospher(wlsCRT >= 390 & wlsCRT <= 780, 2:end);

% calculate xyY coordinates of CRT
xyYCRT = XYZToxyY(T_xyz*rgbCRT);
idxCRT = convhull(xyYCRT(1,:), xyYCRT(2,:));

% calculate photoreceptor activations of CRT
ssCRT = T_cies026*rgbCRT;
idxSSCRT = convhull(ssCRT(1,:), ssCRT(2,:));

%% load in primaries for the Display++
blueDP=load('display++/Blue.mat');
greenDP=load('display++/Green.mat');
redDP=load('display++/Red.mat');
wlsDP = redDP.Lambda;
rgbDP = [redDP.Radiance(wlsDP >=390 & wlsDP <=780)',greenDP.Radiance(wlsDP >=390 & wlsDP <=780)',blueDP.Radiance(wlsDP >=390 & wlsDP <=780)'];
clear blueDP greenDP redDP

% calculate xyY coordinates of Display++
xyYDP = XYZToxyY(T_xyz*rgbDP);
idxDP = convhull(xyYDP(1,:), xyYDP(2,:));

% calculate photoreceptor activations of Display ++
ssDP = T_cies026*rgbDP;
idxSSDP = convhull(ssDP(1,:), ssDP(2,:));

%% load in primaries for the LCD
blueLCD=load('dellLCD/blue.mat');
greenLCD=load('dellLCD/green.mat');
redLCD=load('dellLCD/red.mat');
wlsLCD = redLCD.Lambda;
rgbLCD = [redLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',greenLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',blueLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)'];
clear blueLCD greenLCD redLCD

% calculate xyY coordinates of LCD
xyYLCD = XYZToxyY(T_xyz*rgbLCD);
idxLCD = convhull(xyYLCD(1,:), xyYLCD(2,:));

% calculate photoreceptor activations of LCD
ssLCD = T_cies026*rgbLCD;
idxSSLCD = convhull(ssLCD(1,:), ssLCD(2,:));

%% mythbusting sRGB gamut

% plot on xy diagram
figure('defaultAxesFontSize',18)
subplot(1,2,1)
plotChromaticity();
hold on;
h(1)=plot(xyYCRT(1,idxCRT),xyYCRT(2,idxCRT),'k-','LineWidth',2);
h(2)=plot(xyYDP(1,idxDP),xyYDP(2,idxDP),'b-','LineWidth',2);
h(3) = plot(xyYLCD(1,idxLCD),xyYDP(2,idxLCD),'-','Color',[0.75,0.75,0.75],'LineWidth',2);
legend(h,{'CRT','DP','LCD'});
xlabel('x');
ylabel('y');

% plot primaries out
subplot(1,2,2)
g(1) = plot(390:780, rgbCRT(:,1), 'r-','LineWidth',2);
hold on;
plot(390:780, rgbCRT(:,2), 'g-','LineWidth',2);
plot(390:780, rgbCRT(:,3), 'b-','LineWidth',2);
g(2) = plot(390:780, rgbLCD(:,1), 'r-.','LineWidth',2);
plot(390:780, rgbLCD(:,2), 'g-.','LineWidth',2);
plot(390:780, rgbLCD(:,3), 'b-.','LineWidth',2);
g(3)=plot(390:780, rgbDP(:,1), 'r--','LineWidth',2);
plot(390:780, rgbDP(:,2), 'g--','LineWidth',2);
plot(390:780, rgbDP(:,3), 'b--','LineWidth',2);
legend(g,{'CRT','LCD','DP'});
xlabel('Wavelength (nm)');
ylabel('Radiance (uW/sr/m^2)');
xlim([400,800]);

%%
% figure('defaultAxesFontSize',18)
% subplot(3,1,1)
% l(1) = plot3(ssDP(1,idxSSDP),ssDP(3,idxSSDP),ssDP(5,idxSSDP),'b');
% hold on;
% l(2) = plot3(ssCRT(1,idxSSCRT),ssCRT(3,idxSSCRT),ssCRT(5,idxSSCRT),'k');
% l(3) = plot3(ssLCD(1,idxSSLCD),ssLCD(3,idxSSLCD),ssLCD(5,idxSSLCD),'Color',[0.75,0.75,0.75]);
% xlabel('S');
% ylabel('L');
% zlabel('Mel');
% legend(l,{'DP','CRT','LCD'});
% then plot on MB diagram with melanopsin too!

%% taking xy space into 3 dimensions

% plot for xyY axes
figure('defaultAxesFontSize',18)
subplot(2,2,1)
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),xyYDP(3,idxSSDP),'b');
hold on;
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),xyYCRT(3,idxSSCRT),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),xyYLCD(3,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Y');
legend(l,{'DP','CRT','LCD'});

% plot for xy Mel axes
subplot(2,2,2)
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(5,idxSSDP),'b');
hold on;
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(5,idxSSCRT),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(5,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Mel');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,3)
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(2,idxSSDP)+ssDP(3,idxSSDP),'b');
hold on;
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(2,idxSSDP)+ssCRT(3,idxSSDP),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(2,idxSSDP)+ssLCD(3,idxSSDP),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('L+M');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,4)
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(4,idxSSDP),'b');
hold on;
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(4,idxSSDP),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(4,idxSSDP),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Rod');
legend(l,{'DP','CRT','LCD'});

%% calculating luminance, scotopic luminance, and melanopsin for all points in colour space

c=1;
% calculate all possible rgb combinations for
% crt and lcd and the associated xyY and photorecpetor values
for i=0:0.025:1
    for j =0:0.025:1
        for k=0:0.025:1
            % for CRT
            crtCombo = (i*rgbCRT(:,1))+(j*rgbCRT(:,2))+(k*rgbCRT(:,3));
            xyYCRTcombo(:,c) = XYZToxyY(T_xyz*crtCombo);
            ssCRTcombo(:,c) = T_cies026*crtCombo;
            % for LCD
            lcdCombo = (i*rgbLCD(:,1))+(j*rgbLCD(:,2))+(k*rgbLCD(:,3));
            xyYLCDcombo(:,c) = XYZToxyY(T_xyz*lcdCombo);
            ssLCDcombo(:,c) = T_cies026*lcdCombo;
            % for DP
            dpCombo = (i*rgbDP(:,1))+(j*rgbDP(:,2))+(k*rgbDP(:,3));
            xyYDPcombo(:,c) = XYZToxyY(T_xyz*dpCombo);
            ssDPcombo(:,c) = T_cies026*dpCombo;
            c=c+1;
        end
    end
end

% remove NaNs
xyYCRTcombo(isnan(xyYCRTcombo)) = 0;
xyYLCDcombo(isnan(xyYLCDcombo)) = 0;
xyYCDPcombo(isnan(xyYDPcombo)) = 0;

%% get spectral and daylight locus

slRad = getSpectralLocusSpectra(390:780);
dlRad = getDaylightSpectra;

%% get simulated radiant spectra

simRad = getSimulatedSpectra;

%% set up 5nm spacing colorimetry
wls_xyz = T_xyz(:, 1);
T_xyz_5nm = 683*T_xyz(wls_xyz >= 390 & wls_xyz <= 780, 2:end)';
wls_xyz = wls_xyz(wls_xyz >= 390 & wls_xyz <= 780, 1);
% scale for spds with 5nm spacing
wls_xyz_5nm = wls_xyz(1:5:end);
T_xyz_5nm = T_xyz(:,1:5:end);

% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
% scale for spds with 5nm spacing
wls_cies026_5nm = wlsCIES026(1:5:end);
T_cies026_5nm = T_cies026(:,1:5:end);
% remove Nans
T_cies026_5nm(isnan(T_cies026_5nm)) = 0;

%% calculate xyY and photoreceptor activations of spectral locus

% calculate xyY coordinates of spectral locus
xyYSL = XYZToxyY(T_xyz*slRad);
idxSL = convhull(xyYSL(1,:), xyYSL(2,:));

% calculate photoreceptor activations of spectral locus
ssSL = T_cies026*slRad;

%% calculate xyY and photoreceptor activations of simulated spectra

% calculate xyY coordinates of daylight locus
xyYDL = XYZToxyY(T_xyz_5nm*dlRad);
idxDL = convhull(xyYDL(1,:), xyYDL(2,:));

% calculate photoreceptor activations of daylight locus
ssDL = T_cies026_5nm*dlRad;

%% calculate xyY and photoreceptor activations of simulated spectra

% calculate xyY coordinates of simulated
xyYSim = XYZToxyY(T_xyz_5nm*simRad);

% calculate photoreceptor activations of simulated
ssSim = T_cies026_5nm*simRad;

%% plot heatmaps for photopic, scoptopic luminance, and melanopsin
figure('defaultAxesFontSize',18)

% xyY 
subplot(3,4,1)
scatter(xyYCRTcombo(1,:),xyYCRTcombo(2,:),[],ssCRTcombo(2,:)+ssCRTcombo(3,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'L+M';
xlim([0,1]);
ylim([0,1]);
title('CRT');
xlabel('x');
ylabel('y');


subplot(3,4,2)
scatter(xyYLCDcombo(1,:),xyYLCDcombo(2,:),[],ssLCDcombo(2,:)+ssLCDcombo(3,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'L+M';
xlim([0,1]);
ylim([0,1]);
title('LCD');
xlabel('x');
ylabel('y');


subplot(3,4,3)
scatter(xyYDPcombo(1,:),xyYDPcombo(2,:),[],ssDPcombo(2,:)+ssDPcombo(3,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'L+M';
xlim([0,1]);
ylim([0,1]);
title('Display++');
xlabel('x');
ylabel('y');

subplot(3,4,4)
scatter(xyYSim(1,:),xyYSim(2,:),5,(ssSim(2,:)+ssSim(3,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'L+M';
xlim([0,1]);
ylim([0,1]);
title('Simulated Spectra');
xlabel('x');
ylabel('y');

% xyMel
subplot(3,4,5)
scatter(xyYCRTcombo(1,:),xyYCRTcombo(2,:),[],ssCRTcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

subplot(3,4,6)
scatter(xyYLCDcombo(1,:),xyYLCDcombo(2,:),[],ssLCDcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

subplot(3,4,7)
scatter(xyYDPcombo(1,:),xyYDPcombo(2,:),[],ssDPcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');
zlabel('Y');

subplot(3,4,8)
scatter(xyYSim(1,:),xyYSim(2,:),5,(ssSim(5,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

% xyRod
subplot(3,4,9)
scatter(xyYCRTcombo(1,:),xyYCRTcombo(2,:),[],ssCRTcombo(4,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Rod';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

subplot(3,4,10)
scatter(xyYLCDcombo(1,:),xyYLCDcombo(2,:),[],ssLCDcombo(4,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Rod';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

subplot(3,4,11)
scatter(xyYDPcombo(1,:),xyYDPcombo(2,:),[],ssDPcombo(4,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Rod';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');
zlabel('Y');

subplot(3,4,12)
scatter(xyYSim(1,:),xyYSim(2,:),5,(ssSim(4,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Rod';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');


%% plot difference between lum and mel heatmaps
figure('defaultAxesFontSize',18)

subplot(2,2,1)
scatter(xyYCRTcombo(1,:),xyYCRTcombo(2,:),[],(ssCRTcombo(2,:)+ssCRTcombo(3,:))-ssCRTcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Lum-Mel';
xlim([0,1]);
ylim([0,1]);
title('CRT');
xlabel('x');
ylabel('y');


subplot(2,2,2)
scatter(xyYLCDcombo(1,:),xyYLCDcombo(2,:),[],(ssLCDcombo(2,:)+ssLCDcombo(3,:))-ssLCDcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Lum-Mel';
xlim([0,1]);
ylim([0,1]);
title('LCD');
xlabel('x');
ylabel('y');


subplot(2,2,3)
scatter(xyYDPcombo(1,:),xyYDPcombo(2,:),[],(ssDPcombo(2,:)+ssDPcombo(3,:))-ssDPcombo(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Lum-Mel';
xlim([0,1]);
ylim([0,1]);
title('Display++');
xlabel('x');
ylabel('y');

subplot(2,2,4)
scatter(xyYSim(1,:),xyYSim(2,:),5,(ssSim(2,:)+ssSim(3,:))-ssSim(5,:),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Lum-Mel';
xlim([0,1]);
ylim([0,1]);
xlabel('x');

%% plot as mel/lum 
figure('defaultAxesFontSize',18)

subplot(2,2,1)
scatter(xyYCRTcombo(1,:),xyYCRTcombo(2,:),[],ssCRTcombo(5,:)./(ssCRTcombo(2,:)+ssCRTcombo(3,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel/Lum';
xlim([0,1]);
ylim([0,1]);
title('CRT');
xlabel('x');
ylabel('y');

subplot(2,2,2)
scatter(xyYLCDcombo(1,:),xyYLCDcombo(2,:),[],ssLCDcombo(5,:)./(ssLCDcombo(2,:)+ssLCDcombo(3,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel/Lum';
xlim([0,1]);
ylim([0,1]);
title('LCD');
xlabel('x');
ylabel('y');

subplot(2,2,3)
scatter(xyYDPcombo(1,:),xyYDPcombo(2,:),[],ssDPcombo(5,:)./(ssDPcombo(2,:)+ssDPcombo(3,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel/Lum';
xlim([0,1]);
ylim([0,1]);
title('Display++');
xlabel('x');
ylabel('y');

subplot(2,2,4)
scatter(xyYSim(1,:),xyYSim(2,:),5,ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'filled');
hold on;
plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
c = colorbar;
c.Label.String = 'Mel/Lum';
xlim([0,1]);
ylim([0,1]);
xlabel('x');
ylabel('y');

%% plot on xy diagrams - normalisation question - priamries are in absolute radiance units but illuminants have been normalized

% plot on xy diagram
figure('defaultAxesFontSize',18)
subplot(2,2,1)
plotChromaticity();
hold on;
plot(xyYSim(1,:),xyYSim(2,:),'kx','LineWidth',2);
h(1)=plot(xyYCRT(1,idxCRT),xyYCRT(2,idxCRT),'k-','LineWidth',2);
h(2)=plot(xyYDP(1,idxDP),xyYDP(2,idxDP),'b-','LineWidth',2);
h(3) = plot(xyYLCD(1,idxLCD),xyYDP(2,idxLCD),'-','Color',[0.75,0.75,0.75],'LineWidth',2);
legend(h,{'CRT','DP','LCD'});
xlabel('x');
ylabel('y');

subplot(2,2,2)
scatter3(xyYSim(1,:),xyYSim(2,:),ssSim(2,:)+ssSim(3,:),'kx');
hold on;
plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),ssSL(2,idxSL)+ssSL(3,idxSL),'r-');
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(2,idxSSDP)+ssDP(3,idxSSDP),'b');
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(2,idxSSCRT)+ssCRT(3,idxSSCRT),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(2,idxSSLCD)+ssLCD(3,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('L+M');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,3)
scatter3(xyYSim(1,:),xyYSim(2,:),ssSim(5,:),'kx');
hold on;
plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),ssSL(2,idxSL)+ssSL(3,idxSL),'r-');
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(5,idxSSDP),'b');
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(5,idxSSCRT),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(5,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Mel');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,4)
scatter3(xyYSim(1,:),xyYSim(2,:),ssSim(4,:),'kx');
hold on;
plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),ssSL(2,idxSL)+ssSL(3,idxSL),'r-');
l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(4,idxSSDP),'b');
l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(4,idxSSCRT),'k');
l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(4,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Rod');
legend(l,{'DP','CRT','LCD'});

% difference of simulated and display
% add in spectral locus and daylight locus to all above plots
% move to MacLeod Boynton space with displays

%% correlations

figure('defaultAxesFontSize',18)
subplot(3,2,1)
scatter((ssSim(2,:)+ssSim(3,:)),ssSim(5,:),'k.');
% hold on;
% h(1)=plot(ssCRT(2,idxCRT)+ssCRT(3,idxCRT),ssCRT(5,idxCRT),'r-','LineWidth',2);
% h(2)=plot(ssDP(2,idxCRT)+ssDP(3,idxCRT),ssDP(5,idxCRT),'b-','LineWidth',2);
% h(3)=plot(ssLCD(2,idxCRT)+ssLCD(3,idxCRT),ssLCD(5,idxCRT),'-','Color',[0.85,0.85,0.85],'LineWidth',2);
% legend(h,{'CRT','DP','LCD'});
xlabel('L+M')
ylabel('Mel')

subplot(3,2,2)
scatter((ssSim(2,:)+ssSim(3,:)),ssSim(4,:),'k.');
xlabel('L+M')
ylabel('Rod')

subplot(3,2,3)
scatter(ssSim(4,:),ssSim(5,:),'k.');
xlabel('Rod')
ylabel('Mel')

subplot(3,2,4)
scatter(ssSim(3,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(1,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('L/L+M')
ylabel('S/L+M')

subplot(3,2,5)
scatter(ssSim(3,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('L/L+M')
ylabel('I/L+M')

subplot(3,2,6)
scatter(ssSim(1,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('S/L+M')
ylabel('I/L+M')

%% plot discrimination histograms for L+M vs Mel
figure()
h(1)=histogram(ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'FaceColor',[0,0,0]);
hold on;
h(2)=histogram(ssCRTcombo(5,:)./(ssCRTcombo(2,:)+ssCRTcombo(3,:)),'FaceColor',[1,0,0]);
h(3)=histogram(ssLCDcombo(5,:)./(ssLCDcombo(2,:)+ssLCDcombo(3,:)),'FaceColor',[0.75,0.75,0.75]);
h(4)=histogram(ssDPcombo(5,:)./(ssDPcombo(2,:)+ssDPcombo(3,:)),'FaceColor',[0,0,1]);
legend(h,{'Sim','CRT','LCD','DP'})

%% normalizing display gamuts
%% plot in DKL space with mel axis instead of lum axis and plot display gamuts
%% plot in MB space with extended mel axis and always plot displays