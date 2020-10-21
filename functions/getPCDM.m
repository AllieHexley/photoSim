% function to find the photoreceptor correlation distortions
% created by ACH 01/07/2020
      
function [display] = getPCDM(display,Sim);

% inputs:
% 1) structure of display
% 2) structure of simulated real world spectra

% outputs:
% display structure containing fields for distorted photoreceptor
% correlations

% pairings for correlations
pairs = [1,2;1,3;1,4;1,5;,2,3;2,4;2,5;3,4;3,5;4,5];
pairNames = ['S','M';'S','L';'S','R';'S','I';'M','L';'M','R';'M','I';'L','R';'L','I';'R','I'];

for ii=1:length(pairs)
    [rho{ii}, pval{ii}] = corrcoef(display.ssDistorted(pairs(ii,1),:),display.ssDistorted(pairs(ii,2),:));
   	display.photoreceptorCorrelations(ii) = rho{ii}(1,2);
    display.photoreceptorCorrelationsDistortion(ii) = round(100*((display.photoreceptorCorrelations(ii)-Sim.photoreceptorCorrelations(ii))./Sim.photoreceptorCorrelations(ii)),1);
end

display.correlationLabels = pairNames;

end