%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (390:1:780)';
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
T_cies026 = T_cies026(:,11:end);
% remove Nans
T_cies026(isnan(T_cies026)) = 0;

%% get spectral locus
spectralLocus = getSpectralLocus(T_cies026);

%% get daylight locus
