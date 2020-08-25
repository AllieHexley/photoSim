% function to plot real-world and distorted photoreceptor signals in MacLeod-Boynton space
% created by ACH 01/07/2020
 
clear all;
close all;
clc;
     
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%% plot for each display (CRT,LCD,Display++)

colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5];
labs = ['S';'M';'L';'R';'I'];

figure('defaultAxesFontSize',18)
for i=1:4
    h(1) = plotSimMB(Sim,1);
    h(2) = plotDistortionsMB(CRT,4,colsList,1);
    h(3) = plotDistortionsMB(LCD,7,colsList,2);
    h(4) = plotDistortionsMB(DP,10,colsList,3);
end

function h = plotSimMB(Sim,d)
    subplot(4,3,d)
    lum = (Sim.ss(2,:)+Sim.ss(3,:));
    h=plot(Sim.ss(3,:)./lum,Sim.ss(1,:)./lum,'k.');
    xticklabels({});
    yticklabels({});
    axis square
    ylabel('S/L+M');
    xlabel('L/L+M');

    subplot(4,3,d+1)
    h=plot(Sim.ss(3,:)./lum,Sim.ss(5,:)./lum,'k.');
    xticklabels({});
    yticklabels({});
    axis square
    ylabel('Mel/Lum');
    xlabel('L/L+M');

    subplot(4,3,d+2)
    h=plot(Sim.ss(1,:)./lum,Sim.ss(5,:)./lum,'k.');    xticklabels({});
    xticklabels({});
    yticklabels({});
    axis square
    xlabel('S/L+M');
    ylabel('Mel/Lum');
end

function h = plotDistortionsMB(display,d,colsList,k)
    subplot(4,3,d)
    lum = (display.ssDistorted(2,:)+display.ssDistorted(3,:));
    h=plot(display.ssDistorted(3,:)./lum,display.ssDistorted(1,:)./lum,'.','Color',colsList(k,:));
    xticklabels({});
    yticklabels({});
    axis square
    ylabel('S/L+M');
    xlabel('L/L+M');

    subplot(4,3,d+1)
    h=plot(display.ssDistorted(3,:)./lum,display.ssDistorted(5,:)./lum,'.','Color',colsList(k,:));
    xticklabels({});
    yticklabels({});
    axis square
    ylabel('Mel/Lum');
    xlabel('L/L+M');

    subplot(4,3,d+2)
    h=plot(display.ssDistorted(1,:)./lum,display.ssDistorted(5,:)./lum,'.','Color',colsList(k,:));
    xticklabels({});
    yticklabels({});
    axis square
    xlabel('S/L+M');
    ylabel('Mel/Lum');
end