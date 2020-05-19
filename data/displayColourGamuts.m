% script to convert RGB primaries to xy and xyZ space and to
% MacLeod-Boynton space to see impact of melanopsin on gamut shape
% created by ACH on 13/05/2020

% Set up colorimetry
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

% load in the primaries for a particular display (let's start with the CRT)
load('RGBPhospher.mat');
wlsCRT = RGBPhospher(:,1);
rgbCRT = RGBPhospher(wlsCRT >= 390 & wlsCRT <= 780, 2:end);

% calculate xyY coordinates of CRT
xyYCRT = XYZToxyY(T_xyz*rgbCRT);
idxCRT = convhull(xyYCRT(1,:), xyYCRT(2,:));

% load in primaries for the Display++
blueDP=load('display++/Blue.mat');
greenDP=load('display++/Green.mat');
redDP=load('display++/Red.mat');
wlsDP = redDP.Lambda;
rgbDP = [redDP.Radiance(wlsDP >=390 & wlsDP <=780)',greenDP.Radiance(wlsDP >=390 & wlsDP <=780)',blueDP.Radiance(wlsDP >=390 & wlsDP <=780)'];
clear blueDP greenDP redDP

% calculate xyY coordinates of Display++
xyYDP = XYZToxyY(T_xyz*rgbDP);
idxDP = convhull(xyYDP(1,:), xyYDP(2,:));

% load in primaries for the LCD
blueLCD=load('dellLCD/blue.mat');
greenLCD=load('dellLCD/green.mat');
redLCD=load('dellLCD/red.mat');
wlsLCD = redLCD.Lambda;
rgbLCD = [redLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',greenLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',blueLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)'];
clear blueLCD greenLCD redLCD

% calculate xyY coordinates of Display++
xyYLCD = XYZToxyY(T_xyz*rgbLCD);
idxLCD = convhull(xyYLCD(1,:), xyYLCD(2,:));

% plot on xy diagram
figure('defaultAxesFontSize',18)
plotChromaticity();
hold on;
h(1)=plot(xyYCRT(1,idxCRT),xyYCRT(2,idxCRT),'k-','LineWidth',2);
h(2)=plot(xyYDP(1,idxDP),xyYDP(2,idxDP),'b-','LineWidth',2);
h(3) = plot(xyYLCD(1,idxLCD),xyYDP(2,idxLCD),'-','Color',[0.75,0.75,0.75],'LineWidth',2);
legend(h,{'CRT','DP','LCD'});
xlabel('x');
ylabel('y');

% now plot on MB diagram with just cones
ssDP = T_cies026*rgbDP;
idxSSDP = convhull(ssDP(1,:), ssDP(2,:));
ssCRT = T_cies026*rgbCRT;
idxSSCRT = convhull(ssCRT(1,:), ssCRT(2,:));
ssLCD = T_cies026*rgbLCD;
idxSSLCD = convhull(ssLCD(1,:), ssLCD(2,:));

figure('defaultAxesFontSize',18)
l(1) = plot3(ssDP(1,idxSSDP),ssDP(3,idxSSDP),ssDP(5,idxSSDP),'b');
hold on;
l(2) = plot3(ssCRT(1,idxSSCRT),ssCRT(3,idxSSCRT),ssCRT(5,idxSSCRT),'k');
l(3) = plot3(ssLCD(1,idxSSLCD),ssLCD(3,idxSSLCD),ssLCD(5,idxSSLCD),'Color',[0.75,0.75,0.75]);
xlabel('S');
ylabel('L');
zlabel('Mel');
legend(l,{'DP','CRT','LCD'});
% then plot on MB diagram with melanopsin too!