% script created by ACH 09/04/2020
% script to plot spectral and daylight locus etc
% [S,M,L,R,I]

clear all;
clc;

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

%% get spectral locus
spectralLocus = getSpectralLocus(T_cies026_1nm);

%% get daylight locus
daylightLocus = getDaylightLocus(T_cies026_5nm);

%% get Monte-Carlo simulations of possible illuminants
photoSim = getSimulatedResponses(T_cies026_5nm);

%% plot both for activation space
figure('defaultAxesFontSize',18)
labels = {'S','M','L','R','I'};
count = 1;
whichPlot = 1;
for i=1:5
    for j=1:5
        if ismember(count,[2,3,4,5,8,9,10,14,15,20])==1 
            subplot(3,4,whichPlot)
            plot(photoSim(:,:,i),photoSim(:,:,j),'k.');
            hold on;
            plot(spectralLocus(:,i),spectralLocus(:,j),'r-');
            plot(daylightLocus(:,i),daylightLocus(:,j),'bx');
            set(gca,'YScale','log');
            set(gca,'XScale','log');
            xlabel(labels{i});
            ylabel(labels{j});
            count = count+1;
            whichPlot = whichPlot+1;
        else
            count = count+1;
            continue;
        end
    end
end

%% convert to MacLeod Boynton space
