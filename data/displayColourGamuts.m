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


% display slices here - everything above 0 lum in first plot of slice,
% everything above 10 lum in second slice, etc - which can get from Y value
% also want to do this for melanopsin slices
% LCD up to 291; CRT up to 69; DP up to 132
% for mel: LCD up to 0.3415; CRT up to 0.1006; DP up to 0.1392
% to plot as slice can just set all lum > 10 equal to 10 in plot so forms
% slice in 3d - also plot as full 3d mesh

% define range of luminance levels for slice
lumLevs = [1,10,20,30,40,50,60,100,125,150,200,250];

% define color vector
% set first colour to black
cols = [0, 0, 0];
count = 1./length(lumLevs);
for k=1:length(lumLevs)-1
    cols = [cols; repmat((count),[1,3])];
    count = count + 1./length(lumLevs);
end
% fix colors later

figure('defaultAxesFontSize',18)
sgtitle('Luminance Slices');
for i=1:length(lumLevs)
    % for CRT
    subplot(3,length(lumLevs),i)
    plot(xyYCRTcombo(1,xyYCRTcombo(3,:)>lumLevs(i)),xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('CRT')
    end
    title(num2str(lumLevs(i)));
    % for LCD
    subplot(3,length(lumLevs),length(lumLevs)+i)
    plot(xyYLCDcombo(1,xyYLCDcombo(3,:)>lumLevs(i)),xyYLCDcombo(2,xyYLCDcombo(3,:)>lumLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('LCD');
    end
    % for Display ++
    subplot(3,length(lumLevs),(2*length(lumLevs))+i)
    plot(xyYDPcombo(1,xyYDPcombo(3,:)>lumLevs(i)),xyYDPcombo(2,xyYDPcombo(3,:)>lumLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('DP');
    end
end;

% and for melanopsin
melLevs = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.325];

figure('defaultAxesFontSize',18)
sgtitle('Melanopsin slices');
for i=1:length(melLevs)
    % for CRT
    subplot(3,length(melLevs),i)
    plot(xyYCRTcombo(1,ssCRTcombo(5,:)>melLevs(i)),xyYCRTcombo(2,ssCRTcombo(5,:)>melLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('CRT')
    end
    title(num2str(melLevs(i)));
    % for LCD
    subplot(3,length(melLevs),length(melLevs)+i)
    plot(xyYLCDcombo(1,ssLCDcombo(5,:)>melLevs(i)),xyYLCDcombo(2,ssLCDcombo(5,:)>melLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('LCD');
    end
    % for Display ++
    subplot(3,length(melLevs),(2*length(melLevs))+i)
    plot(xyYDPcombo(1,ssDPcombo(5,:)>melLevs(i)),xyYDPcombo(2,ssDPcombo(5,:)>melLevs(i)),'.','Color',cols(i,:));
    hold on;
    plot(xyYSL(1,idxSL),xyYSL(2,idxSL),'k-');
    xlim([0,1]);
    ylim([0,1]);
    if i==1
        ylabel('DP');
    end
end;

% and for rods

%% and as 3d slices

% i=1;
% % % try with convex hulls
% lowLumCRT = [xyYCRTcombo(1,xyYCRTcombo(3,:)>lumLevs(i));xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i));repmat(lumLevs(i),[1,length(xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i)))])];
% %     
% convCRT = boundary(lowLumCRT(1,:)',lowLumCRT(2,:)');
% figure('defaultAxesFontSize',18)
% plot(xyYCRTcombo(1,convCRT),xyYCRTcombo(2,convCRT));

figure('defaultAxesFontSize',18)

% make 3d slice plots
for i=1:length(lumLevs)
    subplot(2,3,1)
    scatter3(xyYCRTcombo(1,xyYCRTcombo(3,:)>lumLevs(i)),xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,250]);
    title('CRT');
    zlabel('Luminance');
    subplot(2,3,2)
    scatter3(xyYLCDcombo(1,xyYLCDcombo(3,:)>lumLevs(i)),xyYLCDcombo(2,xyYLCDcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYLCDcombo(2,xyYLCDcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,250]);
    title('LCD')
    subplot(2,3,3)
    scatter3(xyYDPcombo(1,xyYDPcombo(3,:)>lumLevs(i)),xyYDPcombo(2,xyYDPcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYDPcombo(2,xyYDPcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,250]);
    title('Display++')
end

% repeat for melanopsin - on same plot for 3d slices
% make 3d slice plots
for i=1:length(melLevs)
    subplot(2,3,4)
    scatter3(xyYCRTcombo(1,ssCRTcombo(5,:)>melLevs(i)),xyYCRTcombo(2,ssCRTcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYCRTcombo(2,ssCRTcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,0.35]);
    zlabel('Mel');
    subplot(2,3,5)
    scatter3(xyYLCDcombo(1,ssLCDcombo(5,:)>melLevs(i)),xyYLCDcombo(2,ssLCDcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYLCDcombo(2,ssLCDcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,0.35]);
    subplot(2,3,6)
    scatter3(xyYDPcombo(1,ssDPcombo(5,:)>melLevs(i)),xyYDPcombo(2,ssDPcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYDPcombo(2,ssDPcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    zlim([0,0.35]);
end

%% and add on simulations

figure('defaultAxesFontSize',18)

% make 3d slice plots
for i=1:length(lumLevs)
    subplot(1,2,1)
    sgtitle('CRT')
    scatter3(xyYCRTcombo(1,xyYCRTcombo(3,:)>lumLevs(i)),xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYCRTcombo(2,xyYCRTcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),xyYSim(3,:),'r.','MarkerSize',0.2);
    zlabel('Luminance');
    zlim([0,70]);
    subplot(1,2,2)
    scatter3(xyYCRTcombo(1,ssCRTcombo(5,:)>melLevs(i)),xyYCRTcombo(2,ssCRTcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYCRTcombo(2,ssCRTcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),ssSim(5,:),'r.','MarkerSize',0.2);
    zlabel('Mel');
    zlim([0,0.12]);
end    

figure('defaultAxesFontSize',18)

% make 3d slice plots
for i=1:length(lumLevs)
    subplot(1,2,1)
    sgtitle('LCD')
    scatter3(xyYLCDcombo(1,xyYLCDcombo(3,:)>lumLevs(i)),xyYLCDcombo(2,xyYLCDcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYLCDcombo(2,xyYLCDcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),xyYSim(3,:),'r.','MarkerSize',0.2);
    zlabel('Luminance');
    zlim([0,300]);
    subplot(1,2,2)
    scatter3(xyYLCDcombo(1,ssLCDcombo(5,:)>melLevs(i)),xyYLCDcombo(2,ssLCDcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYLCDcombo(2,ssLCDcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),ssSim(5,:),'r.','MarkerSize',0.2);
    zlabel('Mel');
    zlim([0,0.35]);
end 

figure('defaultAxesFontSize',18)

% make 3d slice plots
for i=1:length(lumLevs)
    subplot(1,2,1)
    sgtitle('DP')
    scatter3(xyYDPcombo(1,xyYDPcombo(3,:)>lumLevs(i)),xyYDPcombo(2,xyYDPcombo(3,:)>lumLevs(i)),repmat(lumLevs(i),[length(xyYDPcombo(2,xyYDPcombo(3,:)>lumLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),xyYSim(3,:),'r.','MarkerSize',0.2);
    zlabel('Luminance');
    zlim([0,150]);
    subplot(1,2,2)
    scatter3(xyYDPcombo(1,ssDPcombo(5,:)>melLevs(i)),xyYDPcombo(2,ssDPcombo(5,:)>melLevs(i)),repmat(melLevs(i),[length(xyYDPcombo(2,ssDPcombo(5,:)>melLevs(i))),1]),'.','MarkerFaceColor',cols(i,:),'MarkerEdgeColor',cols(i,:));
    hold on;
    plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),repmat(0,[length(xyYSL(1,idxSL)),1]),'k-');
    plot3(xyYSim(1,:),xyYSim(2,:),ssSim(5,:),'r.','MarkerSize',0.2);
    zlabel('Mel');
    zlim([0,0.15]);
end 

% issue here is on basis of scaling

%% old 3D heatmaps

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
% TO DO FIX
%l(1) = plot3(xyYDPcombo(1,idxSSDP),xyYDPcombo(2,idxSSDP),ssDPcombo(2,idxSSDP)+ssDPcombo(3,idxSSDP),'b');
%l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(2,idxSSCRT)+ssCRT(3,idxSSCRT),'k');
%l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(2,idxSSLCD)+ssLCD(3,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('L+M');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,3)
scatter3(xyYSim(1,:),xyYSim(2,:),ssSim(5,:),'kx');
hold on;
plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),ssSL(2,idxSL)+ssSL(3,idxSL),'r-');
%l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(5,idxSSDP),'b');
%l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(5,idxSSCRT),'k');
%l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(5,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Mel');
legend(l,{'DP','CRT','LCD'});

subplot(2,2,4)
scatter3(xyYSim(1,:),xyYSim(2,:),ssSim(4,:),'kx');
hold on;
plot3(xyYSL(1,idxSL),xyYSL(2,idxSL),ssSL(2,idxSL)+ssSL(3,idxSL),'r-');
%l(1) = plot3(xyYDP(1,idxSSDP),xyYDP(2,idxSSDP),ssDP(4,idxSSDP),'b');
%l(2) = plot3(xyYCRT(1,idxSSCRT),xyYCRT(2,idxSSCRT),ssCRT(4,idxSSCRT),'k');
%l(3) = plot3(xyYLCD(1,idxSSLCD),xyYLCD(2,idxSSLCD),ssLCD(4,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('x');
ylabel('y');
zlabel('Rod');
legend(l,{'DP','CRT','LCD'});

% difference of simulated and display
% add in spectral locus and daylight locus to all above plots
% move to MacLeod Boynton space with displays

%% correlations

figure('defaultAxesFontSize',18)
subplot(3,3,1)
scatter((ssSim(2,:)+ssSim(3,:)),ssSim(5,:),'k.');
hold on;
h(1)=plot(ssCRT(2,idxCRT)+ssCRT(3,idxCRT),ssCRT(5,idxCRT),'r-','LineWidth',2);
h(2)=plot(ssDP(2,idxCRT)+ssDP(3,idxCRT),ssDP(5,idxCRT),'b-','LineWidth',2);
h(3)=plot(ssLCD(2,idxCRT)+ssLCD(3,idxCRT),ssLCD(5,idxCRT),'-','Color',[0.85,0.85,0.85],'LineWidth',2);
%legend(h,{'CRT','DP','LCD'});
xlabel('L+M')
ylabel('Mel')

subplot(3,3,2)
scatter((ssSim(2,:)+ssSim(3,:)),ssSim(4,:),'k.');
xlabel('L+M')
ylabel('Rod')

subplot(3,3,3)
scatter(ssSim(4,:),ssSim(5,:),'k.');
xlabel('Rod')
ylabel('Mel')

subplot(3,3,4)
scatter(ssSim(3,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(1,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('L/L+M')
ylabel('S/L+M')

subplot(3,3,5)
scatter(ssSim(3,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('L/L+M')
ylabel('I/L+M')

subplot(3,3,6)
scatter(ssSim(1,:)./(ssSim(2,:)+ssSim(3,:)),ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('S/L+M')
ylabel('I/L+M')

%% plot discrimination histograms for L+M vs Mel

% think of more ways to capture these histograms

figure('defaultAxesFontSize',18)
subplot(4,1,1)
h(1)=histogram(ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'FaceColor',[0,0,0]);
xlim([0,2]);
title('Sim')
xlabel('Mel/L+M');
ylabel('Count');
hold on;
subplot(4,1,2)
h(2)=histogram(ssCRTcombo(5,:)./(ssCRTcombo(2,:)+ssCRTcombo(3,:)),'FaceColor',[1,0,0]);
xlim([0,2]);
title('CRT')
xlabel('Mel/L+M');
ylabel('Count');
subplot(4,1,3)
h(3)=histogram(ssLCDcombo(5,:)./(ssLCDcombo(2,:)+ssLCDcombo(3,:)),'FaceColor',[0.75,0.75,0.75]);
xlim([0,2]);
title('LCD')
xlabel('Mel/L+M');
ylabel('Count');
subplot(4,1,4)
h(4)=histogram(ssDPcombo(5,:)./(ssDPcombo(2,:)+ssDPcombo(3,:)),'FaceColor',[0,0,1]);
xlim([0,2]);
title('DP')
xlabel('Mel/L+M');
ylabel('Count');

%% normalizing display primaries/gamuts - right now just using radiance out of spectrocal
%% 3Dness of gamuts - way of plotting gamuts in 3d is not right - think about how to do this through slices of different Y levels
%% metric could be way of computing spread of histogram - overlap of two histograms

%% plot in DKL space with mel axis instead of lum axis and plot display gamuts

figure('defaultAxesFontSize',18)
subplot(1,2,1)
scatter3((ssSim(3,:)-ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),(ssSim(2,:)+ssSim(3,:)),'k.');
hold on;
% plot3((ssSL(3,:)-ssSL(2,:)),ssSL(1,:)-(ssSL(2,:)+ssSL(3,:)),(ssSL(2,:)+ssSL(3,:)),'r-');
%h(1)=plot3((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),ssCRT(1,idxCRT)-(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),'r-','LineWidth',2);
%h(2)=plot3((ssDP(3,idxCRT)-ssDP(2,idxCRT)),ssDP(1,idxCRT)-(ssDP(2,idxCRT)+ssDP(3,idxCRT)),(ssDP(2,idxCRT)+ssDP(3,idxCRT)),'b-','LineWidth',2);
%h(3)=plot3((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),ssLCD(1,idxCRT)-(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),'-','Color',[0.85,0.85,0.85],'LineWidth',2);
%legend(h,{'CRT','DP','LCD'});
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('L+M');

subplot(1,2,2)
scatter3((ssSim(3,:)-ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),(ssSim(5,:)),'k.');
hold on;
% plot3((ssSL(3,:)-ssSL(2,:)),ssSL(1,:)-(ssSL(2,:)+ssSL(3,:)),(ssSL(2,:)+ssSL(3,:)),'r-');
%h(1)=plot3((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),ssCRT(1,idxCRT)-(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),(ssCRT(5,idxCRT)),'r-','LineWidth',2);
%h(2)=plot3((ssDP(3,idxCRT)-ssDP(2,idxCRT)),ssDP(1,idxCRT)-(ssDP(2,idxCRT)+ssDP(3,idxCRT)),(ssDP(5,idxCRT)),'b-','LineWidth',2);
%h(3)=plot3((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),ssLCD(1,idxCRT)-(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),(ssLCD(5,idxCRT)),'-','Color',[0.85,0.85,0.85],'LineWidth',2);
%legend(h,{'CRT','DP','LCD'});
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('Mel');

%% 2D DKL projections with display gamuts

figure('defaultAxesFontSize',18)
subplot(2,3,1)
scatter((ssSim(3,:)-ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),'k.');
hold on;
k(1)=plot((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),ssCRT(1,idxCRT)-(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),'r-');
k(2)=plot((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),ssLCD(1,idxCRT)-(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),'b-');
k(3)=plot((ssDP(3,idxCRT)-ssDP(2,idxCRT)),ssDP(1,idxCRT)-(ssDP(2,idxCRT)+ssDP(3,idxCRT)),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('L-M')
ylabel('S-(L+M)')

subplot(2,3,2)
scatter((ssSim(3,:)-ssSim(2,:)),(ssSim(3,:)+ssSim(2,:)),'k.');
hold on;
k(1)=plot((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),'r-');
k(2)=plot((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),'b-');
k(3)=plot((ssDP(3,idxCRT)-ssDP(2,idxCRT)),(ssDP(2,idxCRT)+ssDP(3,idxCRT)),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('L-M')
ylabel('L+M')

subplot(2,3,3)
scatter((ssSim(3,:)+ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),'k.');
hold on;
k(1)=plot((ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),ssCRT(1,idxCRT)-(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),'r-');
k(2)=plot((ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),ssLCD(1,idxCRT)-(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),'b-');
k(3)=plot((ssDP(3,idxCRT)+ssDP(2,idxCRT)),ssDP(1,idxCRT)-(ssDP(2,idxCRT)+ssDP(3,idxCRT)),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('L+M')
ylabel('S-(L+M)')

subplot(2,3,4)
scatter((ssSim(3,:)-ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),'k.');
hold on;
k(1)=plot((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),ssCRT(1,idxCRT)-(ssCRT(2,idxCRT)+ssCRT(3,idxCRT)),'r-');
k(2)=plot((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),ssLCD(1,idxCRT)-(ssLCD(2,idxCRT)+ssLCD(3,idxCRT)),'b-');
k(3)=plot((ssDP(3,idxCRT)-ssDP(2,idxCRT)),ssDP(1,idxCRT)-(ssDP(2,idxCRT)+ssDP(3,idxCRT)),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('L-M')
ylabel('S-(L+M)')

subplot(2,3,5)
scatter((ssSim(3,:)-ssSim(2,:)),(ssSim(5,:)),'k.');
hold on;
k(1)=plot((ssCRT(3,idxCRT)-ssCRT(2,idxCRT)),ssCRT(5,idxCRT),'r-');
k(2)=plot((ssLCD(3,idxCRT)-ssLCD(2,idxCRT)),ssLCD(5,idxCRT),'b-');
k(3)=plot((ssDP(3,idxCRT)-ssDP(2,idxCRT)),ssDP(5,idxCRT),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('L-M')
ylabel('Mel')

subplot(2,3,6)
scatter((ssSim(5,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),'k.');
hold on;
k(1)=plot(ssCRT(5,idxCRT),ssCRT(1,idxCRT)-(ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),'r-');
k(2)=plot(ssLCD(5,idxCRT),ssLCD(1,idxCRT)-(ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),'b-');
k(3)=plot(ssDP(5,idxCRT),ssDP(1,idxCRT)-(ssDP(3,idxCRT)+ssDP(2,idxCRT)),'g-');
legend(k,{'CRT','LCD','DP'});
xlabel('Mel')
ylabel('S-(L+M)')

%% plot in DKL space where z axis is mel/lum
figure('defaultAxesFontSize',18)
scatter3((ssSim(3,:)-ssSim(2,:)),ssSim(1,:)-(ssSim(2,:)+ssSim(3,:)),ssSim(5,:)./(ssSim(2,:)+ssSim(3,:)),'k.');
xlabel('L-M')
ylabel('S-(L+M)')
zlabel('I/L+M');

%% set up MB chromaticity

% constants
cl = 0.692839;
cm = 0.349676;
cs = 0.0554786;
% calculate VLambda
vLambda = (cl.*T_cies026(3,:))+(cm.*T_cies026(2,:));
% try to caluclate melanopsin coordinate
imb = T_cies026(5,:)./vLambda;
% calculate melanopsin constant
ci = 1./max(imb);
mb = [T_cies026(1,:).*cs;T_cies026(2,:).*cm;T_cies026(3,:).*cl;T_cies026(5,:).*ci];
mb5 = mb(:,1:5:end);

% calculate mb of spectral locus
mbSL = mb*slRad;
mbSim = mb5*simRad;
mbLCD = mb*rgbLCD;
mbCRT = mb*rgbCRT;
mbDP = mb*rgbDP;
idxSL = convhull(mbSL(1,:), mbSL(2,:));

%% plot in MB space

figure('defaultAxesFontSize',18)
subplot(2,2,1)
scatter(mbSim(3,:)./(mbSim(2,:)+mbSim(3,:)),mbSim(1,:)./(mbSim(2,:)+mbSim(3,:)),'k.');
hold on;
plot(mbSL(3,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),mbSL(1,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),'r-');
h(1)=plot(mbCRT(3,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),mbCRT(1,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),'b-');
h(2)=plot(mbLCD(3,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),mbLCD(1,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),'g-');
h(3)=plot(mbDP(3,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),mbDP(1,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),'c-');
legend(h,{'CRT','LCD','DP'});
ylabel('S/L+M');
xlabel('L/L+M');

subplot(2,2,2)
scatter(mbSim(4,:)./(mbSim(2,:)+mbSim(3,:)),mbSim(1,:)./(mbSim(2,:)+mbSim(3,:)),'k.');
hold on;
plot(mbSL(4,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),mbSL(1,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),'r-');
h(1)=plot(mbCRT(4,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),mbCRT(1,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),'b-');
h(2)=plot(mbLCD(4,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),mbLCD(1,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),'g-');
h(3)=plot(mbDP(4,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),mbDP(1,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),'c-');
legend(h,{'CRT','LCD','DP'});
ylabel('S/L+M');
xlabel('I/L+M');

subplot(2,2,3)
scatter(mbSim(3,:)./(mbSim(2,:)+mbSim(3,:)),mbSim(4,:)./(mbSim(2,:)+mbSim(3,:)),'k.');
hold on;
plot(mbSL(3,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),mbSL(4,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),'r-');
h(1)=plot(mbCRT(3,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),mbCRT(4,idxCRT)./(mbCRT(2,idxCRT)+mbCRT(3,idxCRT)),'b-');
h(2)=plot(mbLCD(3,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),mbLCD(4,idxCRT)./(mbLCD(2,idxCRT)+mbLCD(3,idxCRT)),'g-');
h(3)=plot(mbDP(3,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),mbDP(4,idxCRT)./(mbDP(2,idxCRT)+mbDP(3,idxCRT)),'c-');
legend(h,{'CRT','LCD','DP'});
ylabel('I/L+M');

subplot(2,2,4)
scatter3(mbSim(3,:)./(mbSim(2,:)+mbSim(3,:)),mbSim(1,:)./(mbSim(2,:)+mbSim(3,:)),mbSim(4,:)./(mbSim(2,:)+mbSim(3,:)),'k.');
hold on;
plot3(mbSL(3,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),mbSL(1,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),mbSL(4,idxSL)./(mbSL(2,idxSL)+mbSL(3,idxSL)),'r-');
zlabel('I/L+M');
ylabel('S/L+M');
xlabel('L/L+M');

%% plot in cone activations normalized to peak at unity (MB ideas but not MB normalization)

%% 2D DKL projections with display gamuts

% clean up here to define ssSim(3,:)/xxx as quantity itself

figure('defaultAxesFontSize',18)
subplot(1,3,1)
scatter(ssSim(3,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),ssSim(1,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),'k.');
hold on;
k(1)=plot(ssCRT(3,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),ssCRT(1,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),'r-');
k(2)=plot(ssLCD(3,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),ssLCD(1,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),'g-');
k(3)=plot(ssDP(3,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),ssDP(1,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),'b-');
legend(k,{'CRT','LCD','DP'});
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')

subplot(1,3,2)
scatter(ssSim(3,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),ssSim(5,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),'k.');
hold on;
k(1)=plot(ssCRT(3,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),ssCRT(5,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),'r-');
k(2)=plot(ssLCD(3,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),ssLCD(5,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),'g-');
k(3)=plot(ssDP(3,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),ssDP(5,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),'b-');
legend(k,{'CRT','LCD','DP'});
xlabel('L/L+M+S+R+I')
ylabel('Mel/L+M+S+R+I')

subplot(1,3,3)
scatter(ssSim(5,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),ssSim(1,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),'k.');
hold on;
k(1)=plot(ssCRT(5,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),ssCRT(1,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),'r-');
k(2)=plot(ssLCD(5,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),ssLCD(1,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),'g-');
k(3)=plot(ssDP(5,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),ssDP(1,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),'b-');
legend(k,{'CRT','LCD','DP'});
xlabel('Mel/L+M+S+R+I')
ylabel('S/L+M+S+R+I')


%% remove below plot but try to do a histogram or spread plot to replace this to quantify along each projection, quantify in each space, and quantify in 3D space - quantify volume
figure()
scatter((ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),ssSim(1,:)./(ssSim(5,:)+ssSim(4,:)+ssSim(1,:)+ssSim(3,:)+ssSim(2,:)),'k.');
hold on;
k(1)=plot((ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),ssCRT(1,idxCRT)./(ssCRT(5,idxCRT)+ssCRT(4,idxCRT)+ssCRT(1,idxCRT)+ssCRT(3,idxCRT)+ssCRT(2,idxCRT)),'r-');
k(2)=plot((ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),ssLCD(1,idxCRT)./(ssLCD(5,idxCRT)+ssLCD(4,idxCRT)+ssLCD(1,idxCRT)+ssLCD(3,idxCRT)+ssLCD(2,idxCRT)),'g-');
k(3)=plot((ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),ssDP(1,idxCRT)./(ssDP(5,idxCRT)+ssDP(4,idxCRT)+ssDP(1,idxCRT)+ssDP(3,idxCRT)+ssDP(2,idxCRT)),'b-');
legend(k,{'CRT','LCD','DP'});
xlabel('L+M+S+R+I')
ylabel('S/L+M+S+R+I')
%% try to count number of points in each region - apply this to above plots and other plots that I like
%% figure out how to take all of this into 3D for 3D projections

% try to define triangle polygon from coordinates of CRT gamut
pgonCRT = polyshape(xyYCRT(1,:),xyYCRT(2,:));
figure
plot(pgonCRT);
hold on;
[in, on] = inpolygon(xyYSim(1,:),xyYSim(2,:),xyYCRT(1,:),xyYCRT(2,:));
sum(in)+sum(on)
x=xyYSim(1,:);
y=xyYSim(2,:);
plot(x(in),y(in),'rx');
plot(x(~in),y(~in),'bo');

%% for xyY plots vs sim - make Y and mel L+M/L+M+S+R+I and I/L+M+S+R+I to avoid normalization issues!

%% add in slices for 3D info for the above LM space plots


