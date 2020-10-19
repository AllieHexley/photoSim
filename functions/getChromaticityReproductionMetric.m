% function to find the photoreceptor correlation distortions
% created by ACH 01/07/2020
      
function [display] = getChromaticityReproductionMetric(display,Sim);

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing field for chromaticity reproduction metric

totalSpec = 39699;
[in, on] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),display.xyYMax(1,display.idx),display.xyYMax(2,display.idx));
display.chromaticityReproductionMetric=100.*(sum(in)./totalSpec);

end