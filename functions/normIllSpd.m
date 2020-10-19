% created by ACH 09/04/2020
% script to normalize illuminant spd so integral under illuminant spd = 1;
% as in delta functions which model the spectral locus

function normSpd = normIllSpd(spd, wlsSpd)

for i=1:size(spd,2)
    % calculate integral of illuminant spectra
    A(i) = trapz(wlsSpd, spd(:,i));
    % divide by integral of illuminant spectra so illuminant has unit area
    normSpd(:,i) = spd(:,i)./A(i);
end

end