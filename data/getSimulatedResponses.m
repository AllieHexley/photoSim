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
wlsRef = ref(11:5:401,1);
% resample to match the spectra and cone fundamentals space
ref = ref(11:5:401,2:end);

% normalise illuminant spectra
normSpd = normIllSpd(spd, wlsSpd);

% calculate simulated radiant spectra
k=1;
for i=1:size(normSpd,2)
    for j = 1:size(ref,2)
        rad(:,k) = (ref(:,j).*normSpd(:,i));
        k=k+1;
    end
end

% check norm is correct by summing integrals of k, should equal total
% number of radiance spectra (39699)
normCheck=sum(trapz(wlsSpd, rad(:,:)));
for i=1:size(rad,2)
    if trapz(wlsSpd, rad(:,i)) == 1
        normFails(i) = 0;
    elseif trapz(wlsSpd, rad(:,i)) > 1
        normFails(i) = 1;
    else 
        normFails(i) = -1;
    end
end
% the integrals of all the radiance spectra are below 1

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