% function to plot bar graphs showing metrics for each display
% created by ACH 01/07/2020
 
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%% plot for each display (display,LCD,Display++)

colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5];
labs = ['S';'M';'L';'R';'I'];

figure('defaultAxesFontSize',18)
crtMetrics = getMetrics(CRT);
lcdMetrics = getMetrics(LCD);
dpMetrics = getMetrics(DP);
bar([crtMetrics;lcdMetrics;dpMetrics]);
xticks([1,2,3]);
xticklabels(['   CRT   ';'   LCD   ';'Display++']);
legend({'S Distortion','M Distortion', 'L Distortion','R Distortion','I Distortion', 'Chromaticity Reproduction','Distorted Reproduction','RealWorld Reproduction'});
ylim([0,100]);
ylabel('Metric (%)');
axis square

function displayMetrics = getMetrics(display)

displayMetrics = [display.meanDistortion(1),display.meanDistortion(2),display.meanDistortion(3),display.meanDistortion(4),display.meanDistortion(5),display.chromaticityReproductionMetric, display.distortionReproductionMetric,display.realworldReproductionMetric];

end