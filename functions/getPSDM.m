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
    % calculate distortion for all photoreceptors
    display.distortionMetric(i,:) = 100.*((display.ssDistorted(i,:)-Sim.ss(i,:))./Sim.ss(i,:));
    % calculate mean absolute distortion for each photoreceptor
    display.meanDistortion(i) = mean(abs(display.distortionMetric(i,:)),2);
end

end