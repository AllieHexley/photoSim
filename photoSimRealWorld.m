% function to simulate photoreceptor responses in the real world
% created by ACH 03/03/2020

%% load illuminant spectra dataset from Psychtoolbox
% reference: http://dx.doi.org/10.1364/OE.21.010393

load('spd_houser.mat');

%% load reflectance spectra from the TM-30-15 standard

ref = csvread('99Reflectances.csv');

%% load the photoreceptor spectral sensitivites

[T_cies026, S_cies026] = GetCIES026;

%% resample cone fundamentals to wavelength range of radiance spectra
wlsT = 380:1:780;
wlsSpd = 380:5:780;
for tt = 1:5
    T_photo(tt,:) = interp1(wlsT(:),T_cies026(tt,:),wlsSpd(:));
end

%% calculate the real world photoreceptor responses

%photosim = 