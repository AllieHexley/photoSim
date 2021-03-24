% script to produce the simulate real-world spectra and primaries of a CRT,
% LCD, and Display ++ display, and on two hypothetical five primary displays
% created by ACH 29/06/2020

%% prepare workspace

clear all;
close all;
clc;
addpath('data');
addpath('functions');
if ~exist('plots');
    mkdir('plots');
end
addpath(genpath(pwd));

%% check if file already exists, and ask user if they want to re-run simulation

if exist('photosimReferenceDatabase.mat') == 2
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
    disp('Generating real world reference database and retrieving display calibration data...')
    disp('...this script should only take a few seconds to run...')
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
    
    nb5pR = normpdf(390:780,450,(10./2.355));
    nb5pG = normpdf(390:780,500,(10./2.355));
    nb5pB = normpdf(390:780,550,(10./2.355));
    nb5pC = normpdf(390:780,600,(10./2.355));
    nb5pM = normpdf(390:780,650,(10./2.355));
    wlsnb5p = [390:780];
    rgbcmnb5p = [nb5pR',nb5pG',nb5pB',nb5pC',nb5pM'];
    
    % noramlise so area under primaries is 1
    for i=1:size(rgbcmnb5p,2)
        % calculate integral of illuminant spectra
        A(i) = trapz(wlsnb5p, rgbcmnb5p(:,i));
        rgbcmnb5p(:,i) = rgbcmnb5p(:,i)./A(i);
    end
    
    % calculate xyY coordinates of primaries on max of LCD
    xyYnb5p = XYZToxyY(T_xyz*rgbcmnb5p);
    idxnb5p = convhull(xyYnb5p(1,:), xyYnb5p(2,:));
    
    %% simulate a hypothetical broadband 8-bit 5-primary display
    
    bb5pR = normpdf(390:780,450,(40./2.355));
    bb5pG = normpdf(390:780,500,(40./2.355));
    bb5pB = normpdf(390:780,550,(40./2.355));
    bb5pC = normpdf(390:780,600,(40./2.355));
    bb5pM = normpdf(390:780,650,(40./2.355));
    wlsbb5p = [390:780];
    rgbcmbb5p = [bb5pR',bb5pG',bb5pB',bb5pC',bb5pM'];
    
    % noramlise so area under primaries is 1
    for i=1:size(rgbcmbb5p,2)
        % calculate integral of illuminant spectra
        A(i) = trapz(wlsbb5p, rgbcmbb5p(:,i));
        rgbcmbb5p(:,i) = rgbcmbb5p(:,i)./A(i);
    end
    
    % calculate xyY coordinates of primaries on max of LCD
    xyYbb5p = XYZToxyY(T_xyz*rgbcmbb5p);
    idxbb5p = convhull(xyYbb5p(1,:), xyYbb5p(2,:));
    
    %% get spectral and daylight locus
    
    slRad = getSpectralLocusSpectra(390:780); % get spectral locus from 390:780
    
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
    
    %% set up MacLeod-Boynton chromaticity coordinates
    lScale = 0.69283932; 
    mScale = 0.34967567;
    sScale = 0.05547858;
    % scale factors from CVRL MacLeod & Boynton (1979) 10-deg chromaticity 
    % coordinates based on the Stockman & Sharpe (2000) cone fundamentals: http://www.cvrl.org/    

    mb026(2,:) = T_cies026(2,:)*mScale;
    mb026(3,:) = T_cies026(3,:)*lScale;
    mb026(1,:) = T_cies026(1,:)*sScale;
    
    iScale = 1./(max(T_cies026(5,:)./(mb026(2,:)+mb026(3,:)))); % scale melanopsin spectral sensitivity so that I/L+M peaks at 1
    mb026(5,:) = T_cies026(5,:)*iScale;
    
    % rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
    % scale for spds with 5nm spacing
    mb026_5nm = mb026(:,1:5:end);
    % remove Nans
    mb026_5nm(isnan(mb026_5nm)) = 0;
    mb026(isnan(mb026))=0;
    
    %% calculate xyY and photoreceptor activations of spectral locus
    
    % calculate xyY coordinates of spectral locus
    xyYSL = XYZToxyY(T_xyz*slRad);
    idxSL = convhull(xyYSL(1,:), xyYSL(2,:));
    
    % calculate photoreceptor activations of spectral locus
    ssSL = T_cies026*slRad;
    ssmbSL = mb026*slRad;
    mbSL(1,:) = ssmbSL(1,:)./(ssmbSL(2,:)+ssmbSL(3,:));
    mbSL(2,:) = ssmbSL(3,:)./(ssmbSL(2,:)+ssmbSL(3,:));
    mbSL(3,:) = ssmbSL(5,:)./(ssmbSL(2,:)+ssmbSL(3,:));
    
    %% calculate xyY and photoreceptor activations of simulated spectra
    
    % calculate xyY coordinates of simulated spectra
    xyYSim = XYZToxyY(T_xyz_5nm*simRad);
    
    % calculate photoreceptor activations of simulated spectra
    ssSim = T_cies026_5nm*simRad;
    ssmbSim = mb026_5nm*simRad;
    mbSim(1,:) = ssmbSim(1,:)./(ssmbSim(2,:)+ssmbSim(3,:));
    mbSim(2,:) = ssmbSim(3,:)./(ssmbSim(2,:)+ssmbSim(3,:));
    mbSim(3,:) = ssmbSim(5,:)./(ssmbSim(2,:)+ssmbSim(3,:));
    
    %% calculate photoreceptor correlations of simulated spectra
    pairs = [1,2;1,3;1,4;1,5;,2,3;2,4;2,5;3,4;3,5;4,5];
    pairNames = ['S','M';'S','L';'S','R';'S','I';'M','L';'M','R';'M','I';'L','R';'L','I';'R','I'];
    for i=1:length(pairs)
        [rho{i}, pval{i}] = corrcoef(ssSim(pairs(i,1),:),ssSim(pairs(i,2),:));
        photoreceptorCorrelations(i) = rho{i}(1,2);
    end
    correlationLabels = pairNames;
    
    %% save output
    
    CRT = struct('xyYMax', xyYCRT, 'idx', idxCRT, 'spd', rgbCRT);
    DP = struct('xyYMax', xyYDP, 'idx', idxDP, 'spd', rgbDP);
    LCD = struct('xyYMax', xyYLCD, 'idx', idxLCD, 'spd', rgbLCD);
    nb5p = struct('xyYMax', xyYnb5p, 'idx', idxnb5p, 'spd', rgbcmnb5p);
    bb5p = struct('xyYMax', xyYbb5p, 'idx', idxbb5p, 'spd', rgbcmbb5p);
    Sim = struct('xyY', xyYSim, 'ss', ssSim, 'mb', mbSim, 'photoreceptorCorrelations', photoreceptorCorrelations, 'correlationLabels', correlationLabels);
    SL = struct('xyY', xyYSL, 'idx', idxSL, 'ss', ssSL, 'mb', mbSL);
    
    save('photosimReferenceDatabase.mat','CRT','DP','LCD','nb5p','bb5p','Sim','SL');
    
    % clear messy variables
    clear all;
    % load final struct
    load('photosimReferenceDatabase.mat')
    disp('...done')
%% else load the file
else
    disp('Loading file...')
    load('photosimReferenceDatabase.mat')
    disp('...done')
end
