% created by ACH 09/04/2020
% script to get the spectral locus spectra

function slSPD = getSpectralLocusSpectra(wls)

% simulate spectral locus with delta functions

slSPD = zeros(size(wls,2),size(wls,2));
for ii=1:size(wls,2)
    if ii==1
        slSPD(ii,ii) = 1;
    elseif ii==size(wls,2)
        slSPD(ii,ii) = 1;
    else
    slSPD(ii,ii) = 1;
    end
end;

end

