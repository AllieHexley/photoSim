% script to plot some colour gamut visualizations
% created by ACH on 16/06/2020

clear all;
close all;
clc;

%% load data

load('colorimetry.mat');
load('spectralSensitivities.mat');

%% plot all possible space for simulated spectra
norm = 0;
for i=1:5
    norm = norm+Sim.ss(i,:);
    Sim.ssNorm = Sim.ss./norm;
end

%%

figure('defaultAxesFontSize',18)
subplot(1,3,1)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
hold on;
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')
sgtitle('Unconstrained Simulated Spectra');

subplot(1,3,2)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
hold on;
xlabel('L/L+M+S+R+I')
ylabel('I/L+M+S+R+I')

subplot(1,3,3)
scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
hold on;
xlabel('I/L+M+S+R+I')
ylabel('S/L+M+S+R+I')

%% add in constraint that simulated spectra must be reproducible with display primaries
% i.e. go with constraint that we reproduce cones with priority, and see
% how these changes mel space

% calculate RGB to LMS conversion for the CRT display
% Calculate LMS cone values for each primary in the CRT (i.e. RGB to LMS
% matrix)
crtRGB2LMS = (T_cies026(1:3,:)*CRT.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
crtLMS2RGB = inv(crtRGB2LMS);

% calculate RGB to LMS conversion for the CRT display
% Calculate LMS cone values for each primary in the DP (i.e. RGB to LMS
% matrix)
DPRGB2LMS = (T_cies026(1:3,:)*DP.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
DPLMS2RGB = inv(DPRGB2LMS);

% calculate RGB to LMS conversion for the LCD display
% Calculate LMS cone values for each primary in the LCD (i.e. RGB to LMS
% matrix)
LCDRGB2LMS = (T_cies026(1:3,:)*LCD.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
LCDLMS2RGB = inv(LCDRGB2LMS);
%% for all simulated LMS values, calculate mel value and do quiver plots for how much mel value shifts

for i = 1:length(Sim.ss)
    rgbSettingCRT = crtLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedCRT = (rgbSettingCRT'.*CRT.rgb);
    Sim.ssRealizedCRT(:,i) = T_cies026*(rgbRealizedCRT(:,1)+rgbRealizedCRT(:,2)+rgbRealizedCRT(:,3));
    
    rgbSettingDP = DPLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedDP = (rgbSettingDP'.*DP.rgb);
    Sim.ssRealizedDP(:,i) = T_cies026*(rgbRealizedDP(:,1)+rgbRealizedDP(:,2)+rgbRealizedDP(:,3));
    
    rgbSettingLCD = LCDLMS2RGB*Sim.ss(1:3,i);
    rgbRealizedLCD = (rgbSettingLCD'.*LCD.rgb);
    Sim.ssRealizedLCD(:,i) = T_cies026*(rgbRealizedLCD(:,1)+rgbRealizedLCD(:,2)+rgbRealizedLCD(:,3));
end

normCRT = 0;
normDP = 0;
normLCD = 0;
for i=1:5
    normCRT = normCRT+Sim.ssRealizedCRT(i,:);
    Sim.ssRealizedNormCRT = Sim.ssRealizedCRT./normCRT;
    normDP = normDP+Sim.ssRealizedDP(i,:);
    Sim.ssRealizedNormDP = Sim.ssRealizedDP./normDP;
    normLCD = normLCD+Sim.ssRealizedLCD(i,:);
    Sim.ssRealizedNormLCD = Sim.ssRealizedLCD./normLCD;
end

figure('defaultAxesFontSize',18)
subplot(1,3,1)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
hold on;
scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),'r.');
idxCRT = convhull(CRT.xyY(1,:), CRT.xyY(2,:));
k(1)=plot(CRT.ss(3,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(1,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-');
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')
sgtitle('Constrained Simulated Spectra');

subplot(1,3,2)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
hold on;
scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(5,:),'r.');
k(1)=plot(CRT.ss(3,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(5,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-');
xlabel('L/L+M+S+R+I')
ylabel('I/L+M+S+R+I')
title('CRT')

subplot(1,3,3)
scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
hold on;
scatter(Sim.ssRealizedNormCRT(5,:),Sim.ssRealizedNormCRT(1,:),'r.');
k(1)=plot(CRT.ss(5,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),CRT.ss(1,idxCRT)./(CRT.ss(1,idxCRT)+CRT.ss(2,idxCRT)+CRT.ss(3,idxCRT)+CRT.ss(4,idxCRT)+CRT.ss(5,idxCRT)),'g-');
xlabel('I/L+M+S+R+I')
ylabel('S/L+M+S+R+I')

%% calculate RMS between desired and resulting mel level and plot as histogram of RMS levels

% some error happening here...
a=[Sim.ssRealizedCRT(1,:);Sim.ss(1,:)];
arms = rms(a,1);
b=[Sim.ssRealizedCRT(2,:);Sim.ss(2,:)];
brms = rms(b,1);
c=[Sim.ssRealizedCRT(3,:);Sim.ss(3,:)];
crms = rms(c,1);
d=[Sim.ssRealizedCRT(4,:);Sim.ss(4,:)];
drms = rms(d,1);
e=[Sim.ssRealizedCRT(5,:);Sim.ss(5,:)];
erms = rms(e,1);
figure()
subplot(3,2,1)
histogram(arms);
subplot(3,2,2)
histogram(brms);
subplot(3,2,3)
histogram(crms);
subplot(3,2,4)
histogram(drms);
subplot(3,2,5)
histogram(erms);


% subplot(3,3,4)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormDP(3,:),Sim.ssRealizedNormDP(1,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% sgtitle('Constrained Simulated Spectra');
% 
% subplot(3,3,5)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormDP(3,:),Sim.ssRealizedNormDP(5,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('I/L+M+S+R+I')
% title('DP');
% 
% subplot(3,3,6)
% scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormDP(5,:),Sim.ssRealizedNormDP(1,:),'r.');
% xlabel('I/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% 
% subplot(3,3,7)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormLCD(3,:),Sim.ssRealizedNormLCD(1,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')
% sgtitle('Constrained Simulated Spectra');
% 
% subplot(3,3,8)
% scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormLCD(3,:),Sim.ssRealizedNormLCD(5,:),'r.');
% xlabel('L/L+M+S+R+I')
% ylabel('I/L+M+S+R+I')
% title('LCD');
% 
% subplot(3,3,9)
% scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
% hold on;
% scatter(Sim.ssRealizedNormLCD(5,:),Sim.ssRealizedNormLCD(1,:),'r.');
% xlabel('I/L+M+S+R+I')
% ylabel('S/L+M+S+R+I')

%% plot as quiver plots

figure('defaultAxesFontSize',18)
subplot(1,3,1)
%scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
hold on;
b=100;
quiver(Sim.ssNorm(3,1:b:39600),Sim.ssNorm(1,1:b:39600),Sim.ssRealizedNormCRT(3,1:b:39600)-Sim.ssNorm(3,1:b:39600),Sim.ssRealizedNormCRT(1,1:b:39600)-Sim.ssNorm(1,1:b:39600),0,'MaxHeadSize',0.2);
%scatter(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),'r.','LineWidth',0.5);
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')
sgtitle('Constrained Simulated Spectra');

subplot(1,3,2)
b=800;
quiver(Sim.ssNorm(3,1:b:39600),Sim.ssNorm(5,1:b:39600),Sim.ssRealizedNormCRT(3,1:b:39600)-Sim.ssNorm(3,1:b:39600),Sim.ssRealizedNormCRT(5,1:b:39600)-Sim.ssNorm(5,1:b:39600),0,'MaxHeadSize',0.2);
xlabel('L/L+M+S+R+I')
ylabel('I/L+M+S+R+I')
title('CRT')

subplot(1,3,3)
quiver(Sim.ssNorm(5,1:b:39600),Sim.ssNorm(1,1:b:39600),Sim.ssRealizedNormCRT(5,1:b:39600)-Sim.ssNorm(5,1:b:39600),Sim.ssRealizedNormCRT(1,1:b:39600)-Sim.ssNorm(1,1:b:39600),0,'MaxHeadSize',0.2);
xlabel('I/L+M+S+R+I')
ylabel('S/L+M+S+R+I')

%% plot as 3D surface plot

normCRTCombo = 0;
for i=1:5
    normCRTCombo = normCRTCombo+CRT.ssFullRange(i,:);
    CRT.ssFullRangeNorm = CRT.ssFullRange./normCRTCombo;
end
CRT.ssFullRangeNorm(isnan(CRT.ssFullRangeNorm)==1)=0;

figure('defaultAxesFontSize',18)
scatter3(Sim.ssNorm(3,:),Sim.ssNorm(1,:),Sim.ssNorm(5,:),'k.');
hold on;
scatter3(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),Sim.ssRealizedNormCRT(5,:),'r.');

% attempt to produce mesh
[xsim,ysim] = meshgrid(0:0.01:0.8,0:0.01:0.4);
zsim = griddata(Sim.ssNorm(3,:),Sim.ssNorm(1,:),Sim.ssNorm(5,:),xsim,ysim);
M = mesh(xsim,ysim,zsim);
colormap winter

% attempt to produce mesh for Realized Data
[xcrt,ycrt] = meshgrid(0:0.01:0.8,0:0.01:0.4);
zcrt = griddata(Sim.ssRealizedNormCRT(3,:),Sim.ssRealizedNormCRT(1,:),Sim.ssRealizedNormCRT(5,:),xcrt,ycrt);
mesh(xcrt,ycrt,zcrt);
colormap summer

% attempt to produce mesh for gamut of CRT
[xcrtCombo,ycrtCombo] = meshgrid(0:0.01:0.8,0:0.01:0.4);
zcrtCombo = griddata(CRT.ssFullRangeNorm(3,:),CRT.ssFullRangeNorm(1,:),CRT.ssFullRangeNorm(5,:),xcrtCombo,ycrtCombo);
mesh(xcrtCombo,ycrtCombo,zcrtCombo);

shp = alphaShape(Sim.ssNorm(3,:)',Sim.ssNorm(1,:)',Sim.ssNorm(5,:)');
volume(shp)
shpCRT = alphaShape(Sim.ssRealizedNormCRT(3,:)',Sim.ssRealizedNormCRT(1,:)',Sim.ssRealizedNormCRT(5,:)');
volume(shpCRT)
% this doesn't seem to work the way I want it to...
% perhaps try and look at Rafal's code for producing 3D gamut for MPHDR
% display here
shpCRTCombo = alphaShape(CRT.ssFullRangeNorm(3,:)',CRT.ssFullRangeNorm(1,:)',CRT.ssFullRangeNorm(5,:)');
volume(shpCRTCombo)

% possibly plot on 2D projections in all possible combinations with spread
% along all photoreceptor axes 

%% try with SLI matrix instead of LMS matrix

% calculate RGB to LMS conversion for the CRT display
% Calculate LMS cone values for each primary in the CRT (i.e. RGB to LMS
% matrix)
crtRGB2SLI = (T_cies026([1,3,5],:)*CRT.rgb);

% find inverse to convert from LMS to RGB (i.e. LMS to RGB matrix)
crtLMS2RGBSLI = inv(crtRGB2SLI);


%% do the same thing but with and SLI matrix instead of an LMS matrix

for i = 1:length(Sim.ss)
    rgbSettingCRT = crtLMS2RGBSLI*Sim.ss([1,3,5],i);
    rgbRealizedCRT = (rgbSettingCRT'.*CRT.rgb);
    Sim.ssRealizedSLICRT(:,i) = T_cies026*(rgbRealizedCRT(:,1)+rgbRealizedCRT(:,2)+rgbRealizedCRT(:,3));
end

normCRT = 0;
for i=1:5
    normCRT = normCRT+Sim.ssRealizedSLICRT(i,:);
    Sim.ssRealizedNormSLICRT = Sim.ssRealizedSLICRT./normCRT;
end


figure('defaultAxesFontSize',18)
subplot(1,3,1)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(1,:),'k.');
hold on;
scatter(Sim.ssRealizedNormSLICRT(3,:),Sim.ssRealizedNormSLICRT(1,:),'r.');
xlabel('L/L+M+S+R+I')
ylabel('S/L+M+S+R+I')
sgtitle('Constrained Simulated Spectra');

subplot(1,3,2)
scatter(Sim.ssNorm(3,:),Sim.ssNorm(5,:),'k.');
hold on;
scatter(Sim.ssRealizedNormSLICRT(3,:),Sim.ssRealizedNormSLICRT(5,:),'r.');
xlabel('L/L+M+S+R+I')
ylabel('I/L+M+S+R+I')
title('CRT')

subplot(1,3,3)
scatter(Sim.ssNorm(5,:),Sim.ssNorm(1,:),'k.');
hold on;
scatter(Sim.ssRealizedNormSLICRT(5,:),Sim.ssRealizedNormSLICRT(1,:),'r.');
xlabel('I/L+M+S+R+I')
ylabel('S/L+M+S+R+I')