% script to produce the simulate real-world spectra and primaries of a CRT,
% LCD, and Display ++ display, and on two hypothetical displays
% created by ACH 29/06/2020

%% prepare workspace

clear all;
close all;
clc;
addpath('data');
addpath('functions');

%% check if file already exists, and ask user if they want to re-run simulation

if exist('photosimPhotoreceptorSignals.mat') == 2
    % if file exists, ask user if they just want to load in file
    resimulate = input('Do you want to re-run the simulation (y/n)? ','s');
else
    % if the file doesn't exist, simulate
    resimulate ='y'
end

%% if user wants to resimulate gamuts, or if simulated gamuts file does
% not exist then simulate photoreceptor signals then simulate
% else load the file
    
if resimulate == 'y'

    %% set up colorimetry
    % Get the CIE 2015 10degree XYZ functions
    T_xyz = csvread('data/lin2012xyz10e_1_7sf.csv');
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
    
    %% load in the primaries for the CRT
    load('data/CRT/RGBPhospher.mat');
    wlsCRT = RGBPhospher(:,1);
    rgbCRT = RGBPhospher(wlsCRT >= 390 & wlsCRT <= 780, 2:end);
    
    % calculate xyY coordinates of each individual primary on max for the CRT
    xyYCRT = XYZToxyY(T_xyz*rgbCRT);
    idxCRT = convhull(xyYCRT(1,:), xyYCRT(2,:));
    
    %% load in primaries for the Display++
    blueDP=load('data/Display++/Blue.mat');
    greenDP=load('data/Display++/Green.mat');
    redDP=load('data/Display++/Red.mat');
    wlsDP = redDP.Lambda;
    rgbDP = [redDP.Radiance(wlsDP >=390 & wlsDP <=780)',greenDP.Radiance(wlsDP >=390 & wlsDP <=780)',blueDP.Radiance(wlsDP >=390 & wlsDP <=780)'];
    clear blueDP greenDP redDP
    
    % calculate xyY coordinates of each individual primary on max for the Display++
    xyYDP = XYZToxyY(T_xyz*rgbDP);
    idxDP = convhull(xyYDP(1,:), xyYDP(2,:));
    
    %% load in primaries for the LCD
    blueLCD=load('data/LCD/blue.mat');
    greenLCD=load('data/LCD/green.mat');
    redLCD=load('data/LCD/red.mat');
    wlsLCD = redLCD.Lambda;
    rgbLCD = [redLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',greenLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)',blueLCD.Radiance(wlsLCD >=390 & wlsLCD <=780)'];
    clear blueLCD greenLCD redLCD
    
    % calculate xyY coordinates of each individual primary on max for the  LCD
    xyYLCD = XYZToxyY(T_xyz*rgbLCD);
    idxLCD = convhull(xyYLCD(1,:), xyYLCD(2,:));
    
    %% simulate a hypothetical narrowband 5-primary display
    
    fp1R = normpdf(390:780,450,(10./2.355));
    fp1G = normpdf(390:780,500,(10./2.355));
    fp1B = normpdf(390:780,550,(10./2.355));
    fp1C = normpdf(390:780,600,(10./2.355));
    fp1M = normpdf(390:780,650,(10./2.355));
    wlsFP1 = [390:780];
    rgbcmFP1 = [fp1R',fp1G',fp1B',fp1C',fp1M'];
    
    % noramlise so area under primaries is 1
    for i=1:size(rgbcmFP1,2)
        % calculate integral of illuminant spectra
        A(i) = trapz(wlsFP1, rgbcmFP1(:,i));
        rgbcmFP1(:,i) = rgbcmFP1(:,i)./A(i);
    end
    
    % calculate xyY coordinates of primaries on max of LCD
    xyYFP1 = XYZToxyY(T_xyz*rgbcmFP1);
    idxFP1 = convhull(xyYFP1(1,:), xyYFP1(2,:));
    
    %% simulate a hypothetical broadband 8-bit 5-primary display
    
    fp2R = normpdf(390:780,450,(40./2.355));
    fp2G = normpdf(390:780,500,(40./2.355));
    fp2B = normpdf(390:780,550,(40./2.355));
    fp2C = normpdf(390:780,600,(40./2.355));
    fp2M = normpdf(390:780,650,(40./2.355));
    wlsFP2 = [390:780];
    rgbcmFP2 = [fp2R',fp2G',fp2B',fp2C',fp2M'];
    
    % noramlise so area under primaries is 1
    for i=1:size(rgbcmFP2,2)
        % calculate integral of illuminant spectra
        A(i) = trapz(wlsFP2, rgbcmFP2(:,i));
        rgbcmFP2(:,i) = rgbcmFP2(:,i)./A(i);
    end
    
    % calculate xyY coordinates of primaries on max of LCD
    xyYFP2 = XYZToxyY(T_xyz*rgbcmFP2);
    idxFP2 = convhull(xyYFP2(1,:), xyYFP2(2,:));
    
    %% calculate all possible xyY coordinates and photoreceptor activations from the primaries for 8-bit resolution, assuming linearity and independence of the primaries
    
    c=1;
    d=1;
    bit=255;
    % calculate all possible rgb combinations of displays
    for i=(0:(1./bit):1)
        for j =(0:(1./bit):1)
            for k=(0:(1./bit):1)
                % for CRT
                crtCombo(:,c) = (i*rgbCRT(:,1))+(j*rgbCRT(:,2))+(k*rgbCRT(:,3));
                % for LCD
                lcdCombo(:,c) = (i*rgbLCD(:,1))+(j*rgbLCD(:,2))+(k*rgbLCD(:,3));
                % for DP
                dpCombo(:,c) = (i*rgbDP(:,1))+(j*rgbDP(:,2))+(k*rgbDP(:,3));
                c=c+1;
                % include the below section if you want to simulate the five-primary
                % Warning!: this takes up a lot of memory. I would advise
                % reducing the number of bits when simulating the five primary
                % displays!
                            for l=(0:(1./bit):1)
                                for m=(0:(1./bit):1)
                                % for FP1
                                fp1Combo(:,d) = (i*rgbcmFP1(:,1))+(j*rgbcmFP1(:,2))+(k*rgbcmFP1(:,3))+(l*rgbcmFP1(:,4))+(m*rgbcmFP1(:,5));
                                % for FP2
                                fp2Combo(:,d) = (i*rgbcmFP2(:,1))+(j*rgbcmFP2(:,2))+(k*rgbcmFP2(:,3))+(l*rgbcmFP2(:,4))+(m*rgbcmFP2(:,5));
                                d=d+1;
                                end
                            end
            end
        end
    end
    
    % calculate associated colorimetry
    xyYCRTcombo = XYZToxyY(T_xyz*crtCombo);
    ssCRTcombo = T_cies026*crtCombo;
    xyYLCDcombo = XYZToxyY(T_xyz*lcdCombo);
    ssLCDcombo = T_cies026*lcdCombo;
    xyYDPcombo = XYZToxyY(T_xyz*dpCombo);
    ssDPcombo = T_cies026*dpCombo;
    % Again, uncomment if you want to simulate the five primary displays
    xyYFP1combo = XYZToxyY(T_xyz*fp1Combo);
    ssFP1combo = T_cies026*fp1Combo;
    xyYFP2combo = XYZToxyY(T_xyz*fp2Combo);
    ssFP2combo = T_cies026*fp2Combo;
    
    %% get spectral and daylight locus
    
    slRad = getSpectralLocusSpectra(390:780); % get spectral locus from 390:780
    dlRad = getDaylightSpectra;
    
    %% get simulated radiant spectra
    
    simRad = getSimulatedSpectra;
    
    %% set up 5nm spacing colorimetry (for simulated spectra and daylight locus)
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
    
    %% calculate xyY and photoreceptor activations of daylight locus
    
    % calculate xyY coordinates of daylight locus
    xyYDL = XYZToxyY(T_xyz_5nm*dlRad);
    idxDL = convhull(xyYDL(1,:), xyYDL(2,:));
    
    % calculate photoreceptor activations of daylight locus
    ssDL = T_cies026_5nm*dlRad;
    
    %% calculate xyY and photoreceptor activations of simulated spectra
    
    % calculate xyY coordinates of simulated spectra
    xyYSim = XYZToxyY(T_xyz_5nm*simRad);
    
    % calculate photoreceptor activations of simulated spectra
    ssSim = T_cies026_5nm*simRad;
    
    %% calculate photoreceptor correlations of simulated spectra
    pairs = [1,2;1,3;1,4;1,5;,2,3;2,4;2,5;3,4;3,5;4,5];
    pairNames = ['S','M';'S','L';'S','R';'S','I';'M','L';'M','R';'M','I';'L','R';'L','I';'R','I'];
    for i=1:length(pairs)
        [rho{i}, pval{i}] = corrcoef(Sim.ss(pairs(i,1),:),Sim.ss(pairs(i,2),:));
        Sim.photoreceptorCorrelations(i) = rho{i}(1,2);
    end
    Sim.correlationLabels = pairNames;
    
    %% save output
    
    CRT = struct('xyYMax', xyYCRT, 'idx', idxCRT, 'rgb', rgbCRT, 'xyY', xyYCRTcombo, 'ss', ssCRTcombo);
    DP = struct('xyYMax', xyYDP, 'idx', idxDP, 'rgb', rgbDP, 'xyY', xyYDPcombo, 'ss', ssDPcombo);
    LCD = struct('xyYMax', xyYLCD, 'idx', idxLCD, 'rgb', rgbLCD, 'xyY', xyYLCDcombo, 'ss', ssLCDcombo);
    % include this if you are simulating five-primary displays
    FP1 = struct('xyYMax', xyYFP1, 'idx', idxFP1, 'rgbcm', rgbcmFP1, 'xyY', xyYFP1combo, 'ss', ssFP1combo);
    FP2 = struct('xyYMax', xyYFP2, 'idx', idxFP2, 'rgbcm', rgbcmFP2, 'xyY', xyYFP2combo, 'ss', ssFP2combo);
    Sim = struct('xyY', xyYSim, 'ss', ssSim);
    SL = struct('xyY', xyYSL, 'idx', idxSL, 'ss', ssSL);
    DL = struct('xyY', xyYDL, 'idx', idxDL, 'ss', ssDL);
    
    save('photosimPhotoreceptorSignals.mat','CRT','DP','LCD','Sim','SL','DL','FP1','FP2');
    
%% else load the file
else
    load('photosimPhotoreceptorSignals.mat')
    
end
