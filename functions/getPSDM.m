% function to calculate the photoreceptor distortion metric
% created by ACH 01/07/2020

function [display] = getPSDM(display,Sim)

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor signal
% metric

% loop over all photoreceptors to calculate thier distortion metric
for i=1:5
    % calculate distortion for all photoreceptors if reproducible (i.e.
    % primaries non-zero and non-negative)
    display.distortionMetric(i,display.ssReproducible) = 100.*((display.ssDistorted(i,display.ssReproducible)-Sim.ss(i,display.ssReproducible))./Sim.ss(i,display.ssReproducible));
    % calculate mean of absolute value of distortion for each photoreceptor
    display.meanAbsDistortion(i) = mean(abs(display.distortionMetric(i,display.ssReproducible)),2);
    % calculate mean distortion
    display.meanDistortion(i) = mean((display.distortionMetric(i,display.ssReproducible)),2);
    % calculate standard deviation of distortions
    display.stdDistortion(i) = std((display.distortionMetric(i,display.ssReproducible)),[],2);
end
% calculate overall psdm of display - mean of mean of absolute value of distortions
% for each photoreceptor
display.psdm = mean(display.meanAbsDistortion);

end