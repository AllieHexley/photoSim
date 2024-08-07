% function to find the distorted reproduction metric
% created by ACH 01/07/2020
      
function [display] = getLimitedPSRM(display,Sim);

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor
% correlations

% calculate what percentage of simulated real-world spectra can be
% reproduced on the display
totalSpec = length(Sim.ss);
% check if reproducible (e.g. primaries are non zero and non negative)
ifReproducible = display.ssReproducible;
% check if within 1% of error for each signal
% apply different tolerances for each photoreceptor
% 0.09 for S, 0.02 for L and M, 0.14 for rods, and keep 0.01 for Mel
%withinTolerance = (display.ssDistorted+(display.ssDistorted*0.01)) >= Sim.ss & (display.ssDistorted-(display.ssDistorted*0.01)) <= Sim.ss; % to 1% tolerance
withinTolerance(1,:) = (display.ssDistorted(1,:)+(display.ssDistorted(1,:)*0.02)) >= Sim.ss(1,:) & (display.ssDistorted(1,:)-(display.ssDistorted(1,:)*0.02)) <= Sim.ss(1,:); 
withinTolerance(2,:) = (display.ssDistorted(2,:)+(display.ssDistorted(2,:)*0.02)) >= Sim.ss(2,:) & (display.ssDistorted(2,:)-(display.ssDistorted(2,:)*0.02)) <= Sim.ss(2,:); 
withinTolerance(3,:) = (display.ssDistorted(3,:)+(display.ssDistorted(3,:)*0.09)) >= Sim.ss(3,:) & (display.ssDistorted(3,:)-(display.ssDistorted(3,:)*0.09)) <= Sim.ss(3,:); 
withinTolerance(4,:) = (display.ssDistorted(4,:)+(display.ssDistorted(4,:)*0.14)) >= Sim.ss(4,:) & (display.ssDistorted(4,:)-(display.ssDistorted(4,:)*0.14)) <= Sim.ss(4,:); 
withinTolerance(5,:) = (display.ssDistorted(5,:)+(display.ssDistorted(5,:)*0.01)) >= Sim.ss(5,:) & (display.ssDistorted(5,:)-(display.ssDistorted(5,:)*0.01)) <= Sim.ss(5,:); 

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
