% function to simulate photoreceptor responses in the real world
% created by ACH 03/03/2020

%% load illuminant spectra dataset from Psychtoolbox
% reference: http://dx.doi.org/10.1364/OE.21.010393

spd = csvread('318Illuminants.csv');
wlsSpd = spd(:,1);
spd = spd(:,2:end);

%% load reflectance spectra from the TM-30-15 standard

ref = csvread('99Reflectances.csv');
wlsRef = ref(:,1);
ref = ref(:,2:end);

%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (380:1:780)';
% remove Nans
T_cies026(isnan(T_cies026)) = 0;

%% calculate the real world photoreceptor responses

for i=1:318
    for j=1:99
        for k=1:5
            photosim(i,j,k)=T_cies026(k,:)*(ref(:,j).*spd(:,i));
        end
    end
end

%% remove Nans
photosim(isnan(photosim))=0;

%% calculate MacLeod Boynton responses (all variables are X/L+M)
for i=1:318
    for j=1:99
        L(i,j) = photosim(i,j,3)./((photosim(i,j,3)+photosim(i,j,2)));
        S(i,j) = photosim(i,j,1)./((photosim(i,j,3)+photosim(i,j,2)));
        R(i,j) = photosim(i,j,4)./((photosim(i,j,3)+photosim(i,j,2)));
        I(i,j) = photosim(i,j,5)./((photosim(i,j,3)+photosim(i,j,2)));
    end
end

%% plot MacLeod Boynton space in traditional MacLeod Boynton space
figure()
subplot(2,2,1)
plot(L(:),S(:),'x');
xlabel('L/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,2)
plot(R(:),I(:),'x');
xlabel('R/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,3)
plot(L(:),I(:),'x');
xlabel('L/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,4)
plot(I(:),S(:),'x');
xlabel('Mel/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton in 3D
figure()
scatter3(L(:),S(:),I(:),'x');
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('Mel/L+M');

%% maybe label depending on the reflectance and spd
spdLabels = readtable('318Illuminants_sampleTypes.csv','ReadVariableNames',false);

refLabels = csvread('99Reflectances_sampleTypes.csv');
refLabels = refLabels(1,:); %concatenate because reads in empty rows
%1-nature,2-skin,3-textiles,4-paints,5-plastic,6-printed,7-color system

%% label plots with labels
% there is probably a better way to do this but this will do for now
figure()
for j=1:99
    if refLabels(j)==1
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'kx');
        hold on;
        figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'kx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'kx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'kx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'kx');
        hold on;
    elseif refLabels(j)==2
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'bx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'bx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'bx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'bx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'bx');
        hold on;
    elseif refLabels(j)==3
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'rx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'rx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'rx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'rx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'rx');
        hold on;
    elseif refLabels(j)==4
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'gx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'gx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'gx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'gx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'gx');
        hold on;
    elseif refLabels(j)==5
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'cx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'cx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'cx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'cx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'cx');
        hold on;
    elseif refLabels(j)==6
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'mx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'mx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'mx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'mx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'mx');
        hold on;
    else
        figure(3)
        scatter3(L(:,j),S(:,j),I(:,j),'yx');
        hold on;
                figure(4)
        subplot(2,2,1)
        plot(L(:,j),S(:,j),'yx');
        hold on;
        subplot(2,2,2)
        plot(R(:,j),I(:,j),'yx');
        hold on;
        subplot(2,2,3)
        plot(L(:,j),I(:,j),'yx');
        hold on;
        subplot(2,2,4)
        plot(I(:,j),S(:,j),'yx');
        hold on;
    end
end