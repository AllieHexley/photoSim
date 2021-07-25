% simulate distortions
% created by ACH 06/07/2020

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

%% load in photoreceptor signals
load('suppReferenceDatabase.mat');

%% specify distortion matrix
% i.e. which signals are to be fixed 1=S, 2=M, 3=L, 4=R, 5=I
distortionMatrix = '_ReproduceLMSRI';
matchedSignals = [1,2,3,4,5];

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
    
    % Get the photoreceptor spectral sensitivities - and in the spacing for
    % the Manchester display
    % S, M, L, Rod, Mel
    ss = GetCIES026;
    wlsCIES026_M = (400:1:699)';
    T_cies026_M = ss(:,11:310);
    T_cies026_M(isnan(T_cies026_M)) = 0;
    
    % set up MacLeod-Boynton chromaticity coordinates
    lScale = 0.69283932; 
    mScale = 0.34967567;
    sScale = 0.05547858;
    % scale factors from CVRL MacLeod & Boynton (1979) 10-deg chromaticity 
    % coordinates based on the Stockman & Sharpe (2000) cone fundamentals: http://www.cvrl.org/    

    mb026_M(2,:) = T_cies026_M(2,:)*mScale;
    mb026_M(3,:) = T_cies026_M(3,:)*lScale;
    mb026_M(1,:) = T_cies026_M(1,:)*sScale;
    
    iScale = 1./max(T_cies026_M(5,:)./(mb026_M(2,:)+mb026_M(3,:))); % so I/L+M peaks at 1
    mb026_M(5,:) = T_cies026_M(5,:)*iScale;
    
    % rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
    % scale for spds with 5nm spacing
    mb026_5nm_M = mb026_M(:,1:5:end);
    % remove Nans
    mb026_5nm_M(isnan(mb026_5nm_M)) = 0;
    mb026_M(isnan(mb026_M))=0;
    
    % define smallest bit increment for display
    smallestBit = 1./256;
    
    % get signal distortions for five displays
    disp('Step 1/3: Calculating distorted spectra for each display...');
    disp('...this should take less than a minute...')
    [Man] = getDistortions(1:5,Sim,Man,Man.spd,T_cies026_M, mb026_M, smallestBit,5);
    [MPHDR] = getDistortions(1:5,Sim,MPHDR,MPHDR.spd,T_cies026, mb026, smallestBit,6);
    [Macbook_Pro_2009] = getDistortions(1:3,Sim,Macbook_Pro_2009,Macbook_Pro_2009.spd,T_cies026, mb026, smallestBit,3);
    [Macbook_Pro_2014] = getDistortions(1:3,Sim,Macbook_Pro_2014,Macbook_Pro_2014.spd,T_cies026, mb026, smallestBit,3);
    [Macbook_Air] = getDistortions(1:3,Sim,Macbook_Air,Macbook_Air.spd,T_cies026, mb026, smallestBit,3);
    [Surface_Pro] = getDistortions(1:3,Sim,Surface_Pro,Surface_Pro.spd,T_cies026, mb026, smallestBit,3);
    [NEC] = getDistortions(1:3,Sim,NEC,NEC.spd,T_cies026, mb026, smallestBit,3);
    disp('...done');
    
    % get photoreceptor distortion metrics, PSDM, for five displays
    disp('Step 2/3: Running metrics...');
    disp('...this should take seconds...')
    [Man] = getPSDM(Man,Sim);
    [MPHDR] = getPSDM(MPHDR,Sim);
    [Macbook_Pro_2009] = getPSDM(Macbook_Pro_2009,Sim);
    [Macbook_Pro_2014] = getPSDM(Macbook_Pro_2014,Sim);
    [Macbook_Air] = getPSDM(Macbook_Air,Sim);
    [Surface_Pro] = getPSDM(Surface_Pro,Sim);
    [NEC] = getPSDM(NEC,Sim);
    
    % get photoreceptor correlation distortions
    [Man] = getPCDM(Man,Sim);
    [MPHDR] = getPCDM(MPHDR,Sim);
    [Macbook_Pro_2009] = getPCDM(Macbook_Pro_2009,Sim);
    [Macbook_Pro_2014] = getPCDM(Macbook_Pro_2014,Sim);
    [Macbook_Air] = getPCDM(Macbook_Air,Sim);
    [Surface_Pro] = getPCDM(Surface_Pro,Sim);
    [NEC] = getPCDM(NEC,Sim);
    
    % get real world reproduction metric
    [Man] = getPSRM(Man,Sim);
    [MPHDR] = getPSRM(MPHDR,Sim);
    [Macbook_Pro_2009] = getPSRM(Macbook_Pro_2009,Sim);
    [Macbook_Pro_2014] = getPSRM(Macbook_Pro_2014,Sim);
    [Macbook_Air] = getPSRM(Macbook_Air,Sim);
    [Surface_Pro] = getPSRM(Surface_Pro,Sim);
    [NEC] = getPSRM(NEC,Sim);
    
    % get chromaticity diagram reproduciton metric
    [Man] = getColourGamut(Man,Sim);
    [MPHDR] = getColourGamut(MPHDR,Sim);
    [Macbook_Pro_2009] = getColourGamut(Macbook_Pro_2009,Sim);
    [Macbook_Pro_2014] = getColourGamut(Macbook_Pro_2014,Sim);
    [Macbook_Air] = getColourGamut(Macbook_Air,Sim);
    [Surface_Pro] = getColourGamut(Surface_Pro,Sim);
    [NEC] = getColourGamut(NEC,Sim);
    
    disp('...done');
    %% save output
    disp('Step 3/3: Saving output...');
    disp('...this should take seconds...')
    save(fileName,'Man','MPHDR','Macbook_Pro_2009','Sim','Macbook_Pro_2014','Macbook_Air','Surface_Pro','NEC');
   
    clearvars -except fileName
    
    load(fileName);
    
    clear distortionMatric fileName matchedSignals resimulate
    disp('...done');
else
    % else load the file
    disp('Loading metrics...');
    load(fileName);
    clear distortionMatrix fileName matchedSignals resimulate
    disp('...done');
end