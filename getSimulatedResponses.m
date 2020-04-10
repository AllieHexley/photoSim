% created by ACH 09/04/2020
% script to get the Monte-Carlo simulated spectra

function photosim = getSimulatedResponses(activations);

% reference: http://dx.doi.org/10.1364/OE.21.010393
% load illuminantes
spd = csvread('401Illuminants.csv');
wlsSpd = spd(3:end,1);
% cols are different illuminants, rows are diff wavelengths
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
spd = spd(3:end,2:end);

% load reflectance spectra from the TM-30-15 standard
ref = csvread('99Reflectances.csv');
wlsRef = ref(:,1);
% resample to match the spectra and cone fundamentals space
ref = ref(11:5:401,2:end);

% normalise illuminant spectra
normSpd = normIllSpd(spd, wlsSpd);

% calculate L,M,S,R,I responses
for i=1:size(normSpd,2)
    for j=1:size(ref,2)
        for k=1:size(activations,1)
            photosim(i,j,k)=activations(k,:)*(ref(:,j).*normSpd(:,i));
        end
    end
end
% remove Nans
photosim(isnan(photosim))=0;

end