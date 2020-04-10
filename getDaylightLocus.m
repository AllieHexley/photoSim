% created by ACH 09/04/2020
% script to get the daylight locus

function daylightLocus = getDaylightLocus(activations)

% simulate daylight locus
daylightSPD = csvread('daylightSpectra.csv');
wlsSpd = daylightSPD(3:end,1);
% cols are different illuminants, rows are diff wavelengths
daylightSPD = daylightSPD(3:end,2:end);

% normalise daylight locus
normSpd = normIllSpd(daylightSPD, wlsSpd);

% calculate daylight locus in L,M,S,R,I coordinate space
for m=1:size(daylightSPD,2)
    for k=1:size(activations,1)
        daylightLocus(m,k)=activations(k,:)*(normSpd(:,m));
    end
end
% remove Nans
daylightLocus(isnan(daylightLocus))=0;

end