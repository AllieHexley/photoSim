% created by ACH 09/04/2020
% script to get the spectral locus

function spectralLocus = getSpectralLocus(activations)

% simulate spectral locus with delta functions

slSPD = zeros(size(activations,2),size(activations,2));
for ii=1:size(activations,2)
    slSPD(ii,ii) = 1;
end;

% calculate spectral locus in L,M,S,R,I coordinate space
for m=1:size(activations,2)
    for k=1:size(activations,1)
        spectralLocus(m,k)=activations(k,:)*(slSPD(:,m));
    end
end

% remove Nans
spectralLocus(isnan(spectralLocus))=0;

end

