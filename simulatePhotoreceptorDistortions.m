% simulate distortions
% created by ACH 06/07/2020

%% prepare workspace
clear all;
close all;
clc;
addpath('data');
addpath('functions');

%% load in photoreceptor signals
load('photosimPhotoreceptorSignals.mat');

%% specify distortion matrix
% i.e. which signals are to be fixed - must be single photoreceptor for
% now - I think I can just specify combos in 1st input to get Photoreceptor Signal Distortions function but need to check this
% , but will look to implement different spaces later
distortionMatrix = '_ReproduceLMS';
matchedSignals = [1,2,3];

%% check if the type of distortion already exists
fileName = ['photosimPhotoreceptorDistortions' distortionMatrix, '.mat'];
if exist(fileName)==2
    % if file exists, ask user if they just want to load in file
    resimulate = input('Do you want to re-run the simulation (y/n)? ','s');
else
    % if the file doesn't exist, simulate
    resimulate ='y';
end

%% if user wants to resimulate gamuts, or if simulated gamuts file does
% not exist then simulate photoreceptor signals then simulate
% else load the file
if resimulate =='y'
    
    %% find photoreceptor signal distortions introduced when attempting to preproduce real-world spectra on the display
    
    % Get the photoreceptor spectral sensitivities
    % S, M, L, Rod, Mel
    ss = GetCIES026;
    wlsCIES026 = (390:1:780)';
    T_cies026 = ss(:,11:end);
    T_cies026(isnan(T_cies026)) = 0;
    
    % set up MacLeod-Boynton chromaticity coordinates
    lScale = 0.69283932; 
    mScale = 0.34967567;
    sScale = 0.05547858;
    % scale factors from CVRL MacLeod & Boynton (1979) 10-deg chromaticity 
    % coordinates based on the Stockman & Sharpe (2000) cone fundamentals: http://www.cvrl.org/    

    mb026(2,:) = T_cies026(2,:)*mScale;
    mb026(3,:) = T_cies026(3,:)*lScale;
    mb026(1,:) = T_cies026(1,:)*sScale;
    
    iScale = 1./max(T_cies026(5,:)./(mb026(2,:)+mb026(3,:))); % so I/L+M peaks at 1
    mb026(5,:) = T_cies026(5,:)*iScale;
    
    % rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
    % scale for spds with 5nm spacing
    mb026_5nm = mb026(:,1:5:end);
    % remove Nans
    mb026_5nm(isnan(mb026_5nm)) = 0;
    mb026(isnan(mb026))=0;
    
    % define smallest bit increment for display
    smallestBit = 1./256;
    
    % get signal distortions for five displays
    [CRT] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,CRT,CRT.spd,T_cies026, mb026,smallestBit,3);
    [LCD] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,LCD,LCD.spd,T_cies026, mb026,smallestBit,3);
    [DP] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,DP,DP.spd,T_cies026, mb026, smallestBit,3);
    [FP1] = getPhotoreceptorSignalDistortions(1:5,Sim,FP1,FP1.spd,T_cies026, mb026, smallestBit,5);
    [FP2] = getPhotoreceptorSignalDistortions(1:5,Sim,FP2,FP2.spd,T_cies026, mb026, smallestBit,5);
    
    % get photoreceptor distortion metrics for five displays
    [CRT] = getPhotoreceptorDistortionMetric(CRT,Sim);
    [LCD] = getPhotoreceptorDistortionMetric(LCD,Sim);
    [DP] = getPhotoreceptorDistortionMetric(DP,Sim);
    [FP1] = getPhotoreceptorDistortionMetric(FP1,Sim);
    [FP2] = getPhotoreceptorDistortionMetric(FP2,Sim);
    
    % get photoreceptor distortions across chromaticity space
    [CRT] = getPSDAcrossChromaticities(CRT,Sim);
    [LCD] = getPSDAcrossChromaticities(LCD,Sim);
    [DP] = getPSDAcrossChromaticities(DP,Sim);
    [FP1] = getPSDAcrossChromaticities(FP1,Sim);
    [FP2] = getPSDAcrossChromaticities(FP2,Sim);
    
    % get photoreceptor correlation distortions
    [CRT] = getPhotoreceptorCorrelationDistortions(CRT,Sim);
    [LCD] = getPhotoreceptorCorrelationDistortions(LCD,Sim);
    [DP] = getPhotoreceptorCorrelationDistortions(DP,Sim);
    [FP1] = getPhotoreceptorCorrelationDistortions(FP1,Sim);
    [FP2] = getPhotoreceptorCorrelationDistortions(FP2,Sim);
    
    % get real world reproduction metric
    [CRT] = getRealWorldReproductionMetric(CRT,Sim);
    [LCD] = getRealWorldReproductionMetric(LCD,Sim);
    [DP] = getRealWorldReproductionMetric(DP,Sim);
    [FP1] = getRealWorldReproductionMetric(FP1,Sim);
    [FP2] = getRealWorldReproductionMetric(FP2,Sim);
    
    % get chromaticity diagram reproduciton metric
    [CRT] = getChromaticityReproductionMetric(CRT,Sim);
    [LCD] = getChromaticityReproductionMetric(LCD,Sim);
    [DP] = getChromaticityReproductionMetric(DP,Sim);
    [FP1] = getChromaticityReproductionMetric(FP1,Sim);
    [FP2] = getChromaticityReproductionMetric(FP2,Sim);
    
    %% save output
    
    save(fileName,'CRT','LCD','DP','FP1','FP2','Sim','SL','DL');
   
    clearvars -except fileName
    
    load(fileName);
    
    clear distortionMatric fileName matchedSignals resimulate
else
    % else load the file
    load(fileName);
    clear distortionMatrix fileName matchedSignals resimulate
end