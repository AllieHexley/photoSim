% created by ACH 09/04/2020
% script to get the daylight spectra

function dlNormSpd = getDaylightSpectra

% simulate daylight locus
daylightSPD = csvread('daylightSpectra.csv');
wlsSpd = daylightSPD(3:end,1);
% cols are different illuminants, rows are diff wavelengths
daylightSPD = daylightSPD(3:end,2:end);

% normalise daylight locus
dlNormSpd = normIllSpd(daylightSPD, wlsSpd);

end