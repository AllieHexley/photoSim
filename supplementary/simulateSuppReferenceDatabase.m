% script to produce the simulate real-world spectra and primaries of a CRT,
% LCD, and Display ++ display, and on two hypothetical five primary displays
% created by ACH 29/06/2020

%% prepare workspace

clear all;
close all;
clc;
addpath('data');
addpath('functions');
if ~exist([pwd '/supplementary_plots'],'dir');
    mkdir('supplementary_plots');
end
addpath(genpath(pwd));

%% check if file already exists, and ask user if they want to re-run simulation

if exist('suppReferenceDatabase.mat') == 2
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

    %% and load in the primaries for Oxford MPHDR display
    mphdr = load('MPHDR.mat');
    wls_mphdr = mphdr.wls;
    rgbcmMPHDR = mphdr.mphdrRGBCMY(wls_mphdr >= 390 & wls_mphdr <= 780, 1:end);
    
    % calculate xyY coordinates of each individual primary on max
    xyYMPHDR = XYZToxyY(T_xyz*rgbcmMPHDR);
    idxMPHDR = convhull(xyYMPHDR(1,:), xyYMPHDR(2,:));
    
    %% and for display
    T_xyz = csvread('data/lin2012xyz10e_1_7sf.csv');
    wls_xyz = T_xyz(:, 1);
    T_xyz_M = 683*T_xyz(wls_xyz >= 400 & wls_xyz <= 699, 2:end)';
    wls_xyz_M = wls_xyz(wls_xyz >= 400 & wls_xyz <= 699, 1);
    
    % Get the photoreceptor spectral sensitivities
    % S, M, L, Rod, Mel
    % convert to 400-700nm range!
    ss = GetCIES026;
    wlsCIES026_M = (400:1:699)';
    T_cies026_M = ss(:,11:310);
    T_cies026_M(isnan(T_cies026_M)) = 0;
    
    %% load in the primaries for the Manchester 5 primary display
    manchester = readtable('Manchester_5Primary.xlsx');
    manchester_primaries = table2array(manchester(:,2:6));
    manchester_wls = table2array(manchester(:,1));
    
    wlsMan = manchester_wls(:,1);
    rgbcmMan = manchester_primaries(manchester_wls >= 390 & manchester_wls <= 780, 1:end);
    
    % calculate xyY coordinates of each individual primary on max 
    xyYMan = XYZToxyY(T_xyz_M*rgbcmMan);
    idxMan = convhull(xyYMan(1,:), xyYMan(2,:));
    
    %% get simulated radiant spectra
    
    simRad = getSimulatedSpectra;
    
    %% set up 5nm spacing colorimetry (for simulated spectra and daylight locus)
    T_xyz = csvread('data/lin2012xyz10e_1_7sf.csv');
    wls_xyz = T_xyz(:, 1);
    T_xyz_5nm = 683*T_xyz(wls_xyz >= 390 & wls_xyz <= 780, 2:end)';
    wls_xyz = wls_xyz(wls_xyz >= 390 & wls_xyz <= 780, 1);
    % scale for spds with 5nm spacing
    wls_xyz_5nm = wls_xyz(1:5:end);
    T_xyz_5nm = T_xyz_5nm(:,1:5:end);
    
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
    
    Man = struct('xyYMax', xyYMan, 'idx', idxMan, 'spd', rgbcmMan);
    MPHDR = struct('xyYMax', xyYMPHDR, 'idx', idxMPHDR, 'spd', rgbcmMPHDR);
    Sim = struct('xyY', xyYSim, 'ss', ssSim, 'mb', mbSim, 'photoreceptorCorrelations', photoreceptorCorrelations, 'correlationLabels', correlationLabels);
    
    save('suppReferenceDatabase.mat','Man','Sim','MPHDR');
    
    % clear messy variables
    clear all;
    % load final struct
    load('suppReferenceDatabase.mat')
    disp('...done')
%% else load the file
else
    disp('Loading file...')
    load('suppReferenceDatabase.mat')
    disp('...done')
end
