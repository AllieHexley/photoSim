% function to find the distorted reproduction metric
% created by ACH 01/07/2020
      
function [display] = getPSRM(display,Sim);

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor
% correlations

% calculate what percentage of simulated real-world spectra can be
% reproduced on the display
totalSpec = 39699;
% check if reproducible (e.g. primaries are non zero and non negative)
ifReproducible = display.ssReproducible;
% check if within 1% of error for each signal
withinTolerance = (display.ssDistorted+(display.ssDistorted*0.01)) >= Sim.ss & (display.ssDistorted-(display.ssDistorted*0.01)) <= Sim.ss; % to 1% tolerance
ifWithinTolerance = (sum(withinTolerance(:,:))==5);
% check if within tolerance and reproducible
if size(ifReproducible,1) == size(ifWithinTolerance,1)
    ifMatch = ifWithinTolerance+ifReproducible;
else % catch if the arrays are transposed
    ifMatch = ifWithinTolerance+ifReproducible';
end
%ifMatch = ifWithinTolerance == 1 & sum(display.ssReproducible(:,1)
display.realworldReproductionMetric = 100.*(sum(ifMatch==2)./totalSpec);

end
