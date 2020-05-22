% created by ACH 22/05/2020
% script to plot photoreceptor correlations and correlations of channels

%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (390:1:780)';
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
% scale for spds with 1nm spacing
T_cies026_1nm = T_cies026(:,11:end);
% remove Nans
T_cies026_1nm(isnan(T_cies026_1nm)) = 0;
% scale for spds with 5nm spacing
T_cies026_5nm = T_cies026(:,11:5:end);
% remove Nans
T_cies026_5nm(isnan(T_cies026_5nm)) = 0;
sensitivities = T_cies026_5nm;

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

% calculate L,M,S,R,I responses
for i=1:size(rad,2)
    for k=1:size(sensitivities,1)
       photosim(i,k)=sensitivities(k,:)*(rad(:,i));
    end
end
% remove Nans
photosim(isnan(photosim))=0;

%% reshape to single vector - rework this into function laters
% individual activations
L = photosim(:,3);
S = photosim(:,1);
I = photosim(:,5);
M = photosim(:,2);
R = photosim(:,4);
% photoreceptor responses matrix prm
prm = [S,M,L,R,I];

%%
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',14);
% histograms
subplot(3,2,1)
plot(L./(L+M),S./(L+M),'b.');
xlabel('S/(L+M)');
ylabel('L/(L+M)');

subplot(3,2,2)
plot(L,M,'b.');
xlabel('L');
ylabel('M');

subplot(3,2,3)
plot((L+M),(L-M),'b.');
xlabel('L+M');
ylabel('L-M');

subplot(3,2,4)
plot(L+M,I,'b.');
xlabel('Luminance');
ylabel('Mel');

subplot(3,2,5)
plot(I./(L+M),L+M,'b.');
xlabel('Luminance');
ylabel('Melanopic efficacy of luminance');
