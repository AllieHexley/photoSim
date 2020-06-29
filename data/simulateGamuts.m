% script to produce the simulated spectra and primaries
% created by ACH 29/06/2020

% script to convert RGB primaries to xy and xyZ space and to
% MacLeod-Boynton space to see impact of melanopsin on gamut shape
% created by ACH on 13/05/2020

clear all;
close all;
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

% normalise primaries so full spectra has unit area
A = trapz(390:780, (rgbCRT(:,1)+rgbCRT(:,2)+rgbCRT(:,3)));
for i=1:3
    rgbCRT(:,i) = rgbCRT(:,i)./A;
end

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

% normalise primaries so full spectra has unit area
A = trapz(390:780, (rgbDP(:,1)+rgbDP(:,2)+rgbDP(:,3)));
for i=1:3
    rgbDP(:,i) = rgbDP(:,i)./A;
end

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

% normalise primaries so full spectra has unit area
A = trapz(390:780, (rgbLCD(:,1)+rgbLCD(:,2)+rgbLCD(:,3)));
for i=1:3
    rgbLCD(:,i) = rgbLCD(:,i)./A;
end

% calculate xyY coordinates of LCD
xyYLCD = XYZToxyY(T_xyz*rgbLCD);
idxLCD = convhull(xyYLCD(1,:), xyYLCD(2,:));

% calculate photoreceptor activations of LCD
ssLCD = T_cies026*rgbLCD;
idxSSLCD = convhull(ssLCD(1,:), ssLCD(2,:));

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

%% save 

CRT = struct('xyY', xyYCRT, 'rgb', rgbCRT, 'ss', ssCRT, 'xyYFullRange', xyYCRTcombo, 'ssFullRange', ssCRTcombo);
DP = struct('xyY', xyYDP, 'rgb', rgbDP, 'ss', ssDP, 'xyYFullRange', xyYDPcombo, 'ssFullRange', ssDPcombo);
LCD = struct('xyY', xyYLCD, 'rgb', rgbLCD, 'ss', ssLCD, 'xyYFullRange', xyYLCDcombo, 'ssFullRange', ssLCDcombo);
Sim = struct('xyY', xyYSim, 'ss', ssSim);

save('colorimetryNormPrimaries.mat','CRT','DP','LCD','Sim');