% function to find the photoreceptor signal distortions of real-world
% spectra per chromaticity
% created by ACH 01/07/2020
      
function [display] = getPSDAcrossChromaticities(display,Sim);

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor signals
% of real-world spectra across chromaticity space

% loop over CIExy diagram split into a 100x100 grid and find mean distortion
% for each photoreceptor of all real-world spectra with a chromaticity that
% falls into a particular grid
xStep=0.01;
yStep=0.01;
i=1; j=1;
x=0:0.01:0.99; y = 0:0.01:0.99;
for k =1:5
    for i=1:100
        for j=1:100
            display.meanPhotoreceptorDistortion(k,i,j)=mean(display.distortionMetric(k,Sim.xyY(1,:)>x(i) & Sim.xyY(1,:)<x(i)+xStep & Sim.xyY(2,:)>y(j) & Sim.xyY(2,:)<y(j)+yStep)); 
        end
    end
end
display.meanPhotoreceptorDistortion(isnan(display.meanPhotoreceptorDistortion)==1)=0;
end
