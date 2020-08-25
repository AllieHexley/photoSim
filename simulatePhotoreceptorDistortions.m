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
distortionMatrix = '_ReproduceSLI';
matchedSignals = [1,3,5];

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
    
    % define smallest bit increment for display
    smallestBit = 1./256;
    
    % get signal distortions for five displays
    [CRT] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,CRT,CRT.rgb,T_cies026,smallestBit,3);
    [LCD] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,LCD,LCD.rgb,T_cies026,smallestBit,3);
    [DP] = getPhotoreceptorSignalDistortions(matchedSignals,Sim,DP,DP.rgb,T_cies026,smallestBit,3);
    [FP1] = getPhotoreceptorSignalDistortions(1:5,Sim,FP1,FP1.rgbcm,T_cies026,smallestBit,5);
    [FP2] = getPhotoreceptorSignalDistortions(1:5,Sim,FP2,FP2.rgbcm,T_cies026,smallestBit,5);
    
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
    
    % get distorted photoreceptor signal reproduction metric
    [CRT] = getDistortionReproductionMetric(CRT,Sim);
    [LCD] = getDistortionReproductionMetric(LCD,Sim);
    [DP] = getDistortionReproductionMetric(DP,Sim);
    [FP1] = getDistortionReproductionMetric(FP1,Sim);
    [FP2] = getDistortionReproductionMetric(FP2,Sim);
    
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
    
else
    load(fileName);
end