% created by ACH 26/05/2020
% script to get the Monte-Carlo simulated spectra

function rad = getSimulatedSpectra_MoreReflectances

%% generate daylight spectra
daylightspd4000 = GetDaylightspd(4000,400:700)
daylightspd6500 = GetDaylightspd(6500,400:700)
daylightspd9000 = GetDaylightspd(9000,400:700)
wlsSpd = 400:700;
%% load reflectance spectra from SOCS database and Morimoto dataset

load('reflectance.mat')
load('SOCSDatasets.mat')

% SOCS is 400:10:700, so interpolate
for i=1:length(Reflectance)
    socs_ref(i,:) = SplineSpd([400, 10, 31], Reflectance(i,:)', [400,1,301]);
end

% and just take 400 to 700 range for Morimoto et al
%%
cg_ref = reflectance.Cambridge_Guina(:,11:311);
cu_ref = reflectance.Cambridge_Uganda(:,11:311);
fred_ref = reflectance.Fred(:,11:311);
matsumoto_ref = reflectance.Matsumoto(:,11:311);
brown_ref = reflectance.Brown(:,11:311);
morimoto_ref = reflectance.Morimoto(:,11:311);

ref = [matsumoto_ref; morimoto_ref; brown_ref; fred_ref; cu_ref; cg_ref; socs_ref]'; 

% normalise illuminant spectra
normSpd = [normIllSpd(daylightspd4000, wlsSpd),normIllSpd(daylightspd6500, wlsSpd),normIllSpd(daylightspd9000, wlsSpd)];
normSpd = normSpd(:,[2,4,6]);

% calculate simulated radiant spectra
k=1;
for i=1:size(normSpd,2)
    for j = 1:size(ref,2)
        rad(:,k) = (ref(:,j).*normSpd(:,i));
        k=k+1;
    end
end

end