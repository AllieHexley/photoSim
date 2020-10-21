% simulate distortions
% created by ACH 06/07/2020

%% prepare workspace
clear all;
close all;
clc;
addpath('data');
addpath('functions');

%% load in photoreceptor signals
load('photosimReferenceDatabase.mat');

%% specify distortion matrix
% i.e. which signals are to be fixed - must be single photoreceptor for
% now - I think I can just specify combos in 1st input to get Photoreceptor Signal Distortions function but need to check this
% , but will look to implement different spaces later
distortionMatrix = '_ReproduceLMS';
matchedSignals = [1,2,3];

%% check if the type of distortion already exists
fileName = ['photosimMetrics' distortionMatrix, '.mat'];
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
    [CRT] = getDistortions(matchedSignals,Sim,CRT,CRT.spd,T_cies026, mb026,smallestBit,3);
    [LCD] = getDistortions(matchedSignals,Sim,LCD,LCD.spd,T_cies026, mb026,smallestBit,3);
    [DP] = getDistortions(matchedSignals,Sim,DP,DP.spd,T_cies026, mb026, smallestBit,3);
    [nb5p] = getDistortions(1:5,Sim,nb5p,nb5p.spd,T_cies026, mb026, smallestBit,5);
    [bb5p] = getDistortions(1:5,Sim,bb5p,bb5p.spd,T_cies026, mb026, smallestBit,5);
    
    % get photoreceptor distortion metrics, PSDM, for five displays
    [CRT] = getPSDM(CRT,Sim);
    [LCD] = getPSDM(LCD,Sim);
    [DP] = getPSDM(DP,Sim);
    [nb5p] = getPSDM(nb5p,Sim);
    [bb5p] = getPSDM(bb5p,Sim);
    
    % get photoreceptor correlation distortions
    [CRT] = getPCDM(CRT,Sim);
    [LCD] = getPCDM(LCD,Sim);
    [DP] = getPCDM(DP,Sim);
    [nb5p] = getPCDM(nb5p,Sim);
    [bb5p] = getPCDM(bb5p,Sim);
    
    % get real world reproduction metric
    [CRT] = getPSRM(CRT,Sim);
    [LCD] = getPSRM(LCD,Sim);
    [DP] = getPSRM(DP,Sim);
    [nb5p] = getPSRM(nb5p,Sim);
    [bb5p] = getPSRM(bb5p,Sim);
    
    % get chromaticity diagram reproduciton metric
    [CRT] = getColourGamut(CRT,Sim);
    [LCD] = getColourGamut(LCD,Sim);
    [DP] = getColourGamut(DP,Sim);
    [nb5p] = getColourGamut(nb5p,Sim);
    [bb5p] = getColourGamut(bb5p,Sim);
    
    %% save output
    
    save(fileName,'CRT','LCD','DP','nb5p','bb5p','Sim','SL');
   
    clearvars -except fileName
    
    load(fileName);
    
    clear distortionMatric fileName matchedSignals resimulate
else
    % else load the file
    load(fileName);
    clear distortionMatrix fileName matchedSignals resimulate
end