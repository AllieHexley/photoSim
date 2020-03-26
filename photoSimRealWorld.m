% function to simulate photoreceptor responses in the real world
% created by ACH 03/03/2020

clear all;
clc;

% define parameters
numSpd = 401;
numRef = 99;
numObs = numSpd*numRef;

%% load illuminant spectra dataset from Psychtoolbox
% reference: http://dx.doi.org/10.1364/OE.21.010393

spd = csvread('401Illuminants.csv');
wlsSpd = spd(:,1);
% cols are different illuminants, rows are diff wavelengths
spd = spd(:,2:end);

%% load reflectance spectra from the TM-30-15 standard

ref = csvread('99Reflectances.csv');
wlsRef = ref(:,1);
% resample to match the spectra
ref = ref(1:5:401,2:end);

%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (380:5:780)';
% remove Nans
T_cies026(isnan(T_cies026)) = 0;
% resample to match the spectra
T_cies026 = T_cies026(:,1:5:401);

%% calculate the real world photoreceptor responses

for i=1:numSpd
    for j=1:numRef
        for k=1:5
            photosim(i,j,k)=T_cies026(k,:)*(ref(:,j).*spd(:,i));
        end
    end
end

%% remove Nans
photosim(isnan(photosim))=0;

%% calculate MacLeod Boynton responses (all variables are X/L+M)
for i=1:numSpd
    for j=1:numRef
        Lnorm(i,j) = photosim(i,j,3)./((photosim(i,j,3)+photosim(i,j,2)));
        Mnorm(i,j) = photosim(i,j,2)./((photosim(i,j,3)+photosim(i,j,2)));
        Snorm(i,j) = photosim(i,j,1)./((photosim(i,j,3)+photosim(i,j,2)));
        Rnorm(i,j) = photosim(i,j,4)./((photosim(i,j,3)+photosim(i,j,2)));
        Inorm(i,j) = photosim(i,j,5)./((photosim(i,j,3)+photosim(i,j,2)));
    end
end

%% reshape to single vector - rework this into function laters
% individual activations
L = reshape(photosim(:,:,3),[numObs,1]);
S = reshape(photosim(:,:,1),[numObs,1]);
I = reshape(photosim(:,:,5),[numObs,1]);
M = reshape(photosim(:,:,2),[numObs,1]);
R = reshape(photosim(:,:,4),[numObs,1]);
% photoreceptor responses matrix prm
prm = [S;M;L;R;I];

% Macleod-Boynton activations
Lmb = reshape(Lnorm,[numObs,1]);
Smb = reshape(Snorm,[numObs,1]);
Imb = reshape(Inorm,[numObs,1]);
Mmb = reshape(Mnorm,[numObs,1]);
Rmb = reshape(Rnorm,[numObs,1]);
% macleod-boynton matrix - mbm
mbm = [Smb;Mmb;Lmb;Rmb;Imb];

%% tidy up a bit
clear i j k 

%% add in spectral locus
% this will be for EEW so need to do this for each individual wavelength to
% get full locus
% set up spectral locus vector
slSPD = zeros(81,81);
for ii=1:81
    slSPD(ii,ii) = 1;
end;
% calculate spectral locus in L,M,S,R,I coordinate space
% remove Nans
spectralLocus(isnan(spectralLocus))=0;
for m=1:81
    for k=1:5
        spectralLocus(m,k)=T_cies026(k,:)*(slSPD(:,m));
    end
end

%% calculate macleod boynton of spectral locus
% for all need to shift so starts at 390 or else will reach infintiy
for n=1:81
    Lsl(n) = spectralLocus(n,3)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Msl(n) = spectralLocus(n,2)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Ssl(n) = spectralLocus(n,1)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Rsl(n) = spectralLocus(n,4)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Isl(n) = spectralLocus(n,5)./((spectralLocus(n,3)+spectralLocus(n,2)));
end
Lsl(isnan(Lsl))=0;
Msl(isnan(Msl))=0;
Ssl(isnan(Ssl))=0;
Ssl = Ssl./(max(Ssl));
Rsl(isnan(Rsl))=0;
Isl(isnan(Isl))=0;
%% plot spectral locus alone
figure()
subplot(2,2,1)
plot(Lsl(:),Ssl(:),'ro');
xlabel('L/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,2)
plot(Rsl(:),Isl(:),'ro');
xlabel('R/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,3)
plot(Lsl(:),Isl(:),'ro');
xlabel('L/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,4)
plot(Isl(:),Ssl(:),'ro');
xlabel('Mel/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton in 3D
figure()
scatter3(Lsl(:),Ssl(:),Isl(:),'ro');
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('Mel/L+M');

%% plot spectral locus with real world spectra
% try to plot range of illuminants with the perfect reflectance spectra
figure()
subplot(2,2,1)
plot(Lmb(:),Smb(:),'x');
hold on;
plot(Lsl(:),Ssl(:),'ro');
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


%% do daylight locus assuming perfect reflectance spectra to start

%% then do daylight with actual reflectance spectra

%% then do all illuminants with perfect reflectance spectra

%% then do all illuminants with actual reflectance spectra

%% plot MacLeod Boynton space in traditional MacLeod Boynton space
plotAllPhotoSim(Lmb,Smb,Imb,Rmb);

%% add in daylight locus

%% plot individual activation correlations
figure()
plot(I,I,'x');

%% label depending on the illuminant spectra
[spdLabels, spdLabelsValues] = labelSpd;
    
%% label depending on the reflectance spectra
[refLabels, refLabelsKeys] = labelRef;

%% find correlations
figure()
subplot(3,2,1)
plot(L(:),S(:),'x');
xlabel('L/L+M');
ylabel('S/L+M');
[r,p]=corrcoef([Lcf,Icf]);
text(0.65,0.9,['R=' num2str(r(2))]);
text(0.65,0.8,['p=' num2str(p(2))]);

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(3,2,2)
plot(R(:),I(:),'x');
xlabel('R/L+M');
ylabel('Mel/L+M');
[r,p]=corrcoef([Rcf,Icf]);
text(0.2,0.9,['R=' num2str(r(2))]);
text(0.2,0.8,['p=' num2str(p(2))]);

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(3,2,3)
plot(L(:),I(:),'x');
xlabel('L/L+M');
ylabel('Mel/L+M');
[r,p]=corrcoef([Lcf,Icf]);
text(0.65,0.9,['R=' num2str(r(2))]);
text(0.65,0.8,['p=' num2str(p(2))]);

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(3,2,4)
plot(I(:),S(:),'x');
xlabel('Mel/L+M');
ylabel('S/L+M');
[r,p]=corrcoef([Icf,Scf]);
text(0.2,0.9,['R=' num2str(r(2))]);
text(0.2,0.8,['p=' num2str(p(2))]);

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(3,2,5)
plot(L(:),M(:),'x');
xlabel('L/L+M');
ylabel('M/L+M');
[r,p]=corrcoef([Lcf,Mcf]);
text(0.65,0.55,['R=' num2str(r(2))]);
text(0.65,0.5,['p=' num2str(p(2))]);

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(3,2,6)
plot(L(:),R(:),'x');
xlabel('L/L+M');
ylabel('Rod/L+M');
[r,p]=corrcoef([Lcf,Rcf]);
text(0.65,0.9,['R=' num2str(r(2))]);
text(0.65,0.8,['p=' num2str(p(2))]);

%% plot histograms of photoreceptor responses in natural world not normalised
figure()
sgtitle('Histograms of Activations')
subplot(3,2,1)
histogram(photosim(:,:,3),'FaceColor','r');
title('L');
xlim([0,80]);
ylim([0,3500]);
subplot(3,2,2)
histogram(photosim(:,:,2),'FaceColor','g');
title('M');
xlim([0,80]);
ylim([0,3500]);
subplot(3,2,3)
histogram(photosim(:,:,1),'FaceColor','b');
title('S');
xlim([0,80]);
ylim([0,3500]);
subplot(3,2,4)
histogram(photosim(:,:,4),'FaceColor','m');
title('Rod');
xlim([0,80]);
ylim([0,3500]);
subplot(3,2,5)
histogram(photosim(:,:,5),'FaceColor','c');
title('Mel');
xlim([0,80]);
ylim([0,3500]);

%% plot histograms of photoreceptor responses in natural world
figure()
sgtitle('Histograms of Activations')
subplot(3,2,1)
histogram(L(:),'FaceColor','r');
title('L/L+M');
xlim([0,1]);
ylim([0,1800]);
subplot(3,2,2)
histogram(M(:),'FaceColor','g');
title('M/L+M');
xlim([0,1]);
ylim([0,1800]);
subplot(3,2,3)
histogram(S(:),'FaceColor','b');
title('S/L+M');
xlim([0,1]);
ylim([0,1800]);
subplot(3,2,4)
histogram(R(:),'FaceColor','m');
title('Rod/L+M');
xlim([0,1]);
ylim([0,1800]);
subplot(3,2,5)
histogram(I(:),'FaceColor','c');
title('Mel/L+M');
xlim([0,1]);
ylim([0,1800]);

%%
figure()
sgtitle('Histograms of Activations')
h= histogram(L(:),'FaceColor','r');
hold on;
histogram(M(:),'FaceColor','g');
histogram(S(:),'FaceColor','b');
histogram(R(:),'FaceColor','m');
histogram(I(:),'FaceColor','c');

%% histogram of daylight mel spectra
figure()
subplot(3,2,1)
h1=histogram(L(spdLabelsKeys<10,:));
hold on;
h2=histogram(L(spdLabelsKeys>9,:))
legend([h1,h2],{'Other illuminant','Daylight illuminant'});
title('L/L+M');

subplot(3,2,2)
h1=histogram(M(spdLabelsKeys<10,:));
hold on;
h2=histogram(M(spdLabelsKeys>9,:))
legend([h1,h2],{'Other illuminant','Daylight illuminant'});
title('M/L+M');

subplot(3,2,3)
h1=histogram(S(spdLabelsKeys<10,:));
hold on;
h2=histogram(S(spdLabelsKeys>9,:))
legend([h1,h2],{'Other illuminant','Daylight illuminant'});
title('S/L+M');

subplot(3,2,4)
h1=histogram(R(spdLabelsKeys<10,:));
hold on;
h2=histogram(R(spdLabelsKeys>9,:))
legend([h1,h2],{'Other illuminant','Daylight illuminant'});
title('R/L+M');

subplot(3,2,5)
h1=histogram(I(spdLabelsKeys<10,:));
hold on;
h2=histogram(I(spdLabelsKeys>9,:))
legend([h1,h2],{'Other illuminant','Daylight illuminant'});
title('I/L+M');
%% reshape for curve fitting
Lcf = reshape(L,[numObs,1]);
Scf = reshape(S,[numObs,1]);
Icf = reshape(I,[numObs,1]);

figure()
subplot(2,2,1)
[mbFit1,gof1] = fit([Lcf,Scf],Icf,'poly11');
plot(mbFit1,[Lcf,Scf],Icf);
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');
text(0.4,0.1,['R^2: ', num2str(gof1.rsquare)]);
title('1st order polynomial');

subplot(2,2,2)
[mbFit2,gof2] = fit([Lcf,Scf],Icf,'poly22');
plot(mbFit2,[Lcf,Scf],Icf);
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');
text(0.4,0.1,['R^2: ', num2str(gof2.rsquare)]);
title('2nd order polynomial');

subplot(2,2,3)
[mbFit3,gof3] = fit([Lcf,Scf],Icf,'poly33');
plot(mbFit3,[Lcf,Scf],Icf);
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');
text(0.4,0.1,['R^2: ', num2str(gof3.rsquare)]);
title('3rd order polynomial');

%% attempt PCA
numObs = 318*99;
figure()
X = [reshape(L,[numObs,1]),reshape(S,[numObs,1]),reshape(I,[numObs,1])];
[coeff,score,latent,tsquared,explained] = pca(X);
scatter(score(:,1),score(:,2));
xlabel('1st Principal Component')
ylabel('2nd Principal Component')
text(-0.2,0.4,['1st PC:', num2str(explained(1)),'--' '2nd PC: ', num2str(explained(2))]);


%% label plots with labels
% plot depending on the illuminant spectra label
figure()
cols = [0.1,0.1,0.1;0.2,0.2,0.2;0.3,0.3,0.3;0.4,0.4,0.4;0.5,0.5,0.5;0.6,0.6,0.6;0.7,0.7,0.7;0.8,0.8,0.8;0.9,0.9,0.9;1,0.5,0;1,0.1,0.1;0.1,1,0.1;0.1,0.1,1;0.1,1,1;1,0.1,1;1,1,0.1];
for m=1:16
    subplot(2,1,1)
    plot(L((spdLabelsKeys==m),:),I((spdLabelsKeys==m),:),'x','Color',cols(m,:));
    hold on;
end
subplot(2,1,2)
plot(L(:),I(:),'kx')

%%

% label daylight series
figure()
ii=3
subplot(3,2,1)
plot(L(:,1),S(:,1),'k.');
hold on;
daylightColours = [0.2,1,1;0.2,0.8,1;0.2,0.6,1.0;0.2,0.4,1.0;0.2,0.2,1.0;0.2,0.0,1.0;0.2,0.0,0.8];
p1=plot(L((spdLabelsKeys==10),ii),S((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),ii),S((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),ii),S((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),ii),S((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),ii),S((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),ii),S((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),ii),S((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('L/L+M');
ylabel('S/L+M');
sgtitle('Skin Reflectance Spectra Example');

subplot(3,2,2)
plot(R(:,1),I(:,1),'k.');
hold on;
p1=plot(R((spdLabelsKeys==10),ii),I((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(R((spdLabelsKeys==11),ii),I((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(R((spdLabelsKeys==12),ii),I((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(R((spdLabelsKeys==13),ii),I((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(R((spdLabelsKeys==14),ii),I((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(R((spdLabelsKeys==15),ii),I((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(R((spdLabelsKeys==16),ii),I((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('R/L+M');
ylabel('Mel/L+M');
%title('Natural Reflectance Spectra');

subplot(3,2,3)
plot(L(:,1),I(:,1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),ii),I((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),ii),I((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),ii),I((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),ii),I((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),ii),I((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),ii),I((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),ii),I((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('L/L+M');
ylabel('Mel/L+M');
%title('Natural Reflectance Spectra');

subplot(3,2,4)
plot(S(:,1),I(:,1),'k.');
hold on;
p1=plot(S((spdLabelsKeys==10),ii),I((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(S((spdLabelsKeys==11),ii),I((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(S((spdLabelsKeys==12),ii),I((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(S((spdLabelsKeys==13),ii),I((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(S((spdLabelsKeys==14),ii),I((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(S((spdLabelsKeys==15),ii),I((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(S((spdLabelsKeys==16),ii),I((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('S/L+M');
ylabel('Mel/L+M');
%title('Natural Reflectance Spectra');

subplot(3,2,5)
plot(L(:,1),M(:,1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),ii),M((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),ii),M((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),ii),M((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),ii),M((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),ii),M((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),ii),M((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),ii),M((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('L/L+M');
ylabel('M/L+M');
%title('Natural Reflectance Spectra');

subplot(3,2,6)
plot(L(:,1),R(:,1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),ii),R((spdLabelsKeys==10),ii),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),ii),R((spdLabelsKeys==11),ii),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),ii),R((spdLabelsKeys==12),ii),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),ii),R((spdLabelsKeys==13),ii),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),ii),R((spdLabelsKeys==14),ii),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),ii),R((spdLabelsKeys==15),ii),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),ii),R((spdLabelsKeys==16),ii),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','6500K','7000K','7500K'});
xlabel('L/L+M');
ylabel('R/L+M');
%title('Natural Reflectance Spectra');

%% plots of all daylight illuminants for nature vs other spectra
figure()
subplot(3,2,1)
plot(L((spdLabelsKeys>9),refLabels~=1),S((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
daylightColours = [0.2,1,1;0.2,0.8,1;0.2,0.6,1.0;0.2,0.4,1.0;0.2,0.2,1.0;0.2,0.0,1.0;0.2,0.0,0.8];
p1=plot(L((spdLabelsKeys==10),refLabels==1),S((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),refLabels==1),S((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),refLabels==1),S((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),refLabels==1),S((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),refLabels==1),S((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),refLabels==1),S((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),refLabels==1),S((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('L/L+M');
ylabel('S/L+M');
sgtitle('Daylight Illuminant Spectra Example');

subplot(3,2,2)
plot(R((spdLabelsKeys>9),refLabels~=1),I((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
p1=plot(R((spdLabelsKeys==10),refLabels==1),I((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(R((spdLabelsKeys==11),refLabels==1),I((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(R((spdLabelsKeys==12),refLabels==1),I((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(R((spdLabelsKeys==13),refLabels==1),I((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(R((spdLabelsKeys==14),refLabels==1),I((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(R((spdLabelsKeys==15),refLabels==1),I((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(R((spdLabelsKeys==16),refLabels==1),I((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('R/L+M');
ylabel('I/L+M');

subplot(3,2,3)
plot(L((spdLabelsKeys>9),refLabels~=1),I((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),refLabels==1),I((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),refLabels==1),I((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),refLabels==1),I((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),refLabels==1),I((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),refLabels==1),I((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),refLabels==1),I((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),refLabels==1),I((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('L/L+M');
ylabel('I/L+M');

subplot(3,2,4)
plot(S((spdLabelsKeys>9),refLabels~=1),I((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
p1=plot(S((spdLabelsKeys==10),refLabels==1),I((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(S((spdLabelsKeys==11),refLabels==1),I((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(S((spdLabelsKeys==12),refLabels==1),I((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(S((spdLabelsKeys==13),refLabels==1),I((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(S((spdLabelsKeys==14),refLabels==1),I((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(S((spdLabelsKeys==15),refLabels==1),I((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(S((spdLabelsKeys==16),refLabels==1),I((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('S/L+M');
ylabel('I/L+M');

subplot(3,2,5)
plot(L((spdLabelsKeys>9),refLabels~=1),M((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),refLabels==1),M((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),refLabels==1),M((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),refLabels==1),M((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),refLabels==1),M((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),refLabels==1),M((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),refLabels==1),M((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),refLabels==1),M((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('L/L+M');
ylabel('M/L+M');

subplot(3,2,6)
plot(L((spdLabelsKeys>9),refLabels~=1),R((spdLabelsKeys>9),refLabels~=1),'k.');
hold on;
p1=plot(L((spdLabelsKeys==10),refLabels==1),R((spdLabelsKeys==10),refLabels==1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(L((spdLabelsKeys==11),refLabels==1),R((spdLabelsKeys==11),refLabels==1),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(L((spdLabelsKeys==12),refLabels==1),R((spdLabelsKeys==12),refLabels==1),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(L((spdLabelsKeys==13),refLabels==1),R((spdLabelsKeys==13),refLabels==1),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(L((spdLabelsKeys==14),refLabels==1),R((spdLabelsKeys==14),refLabels==1),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(L((spdLabelsKeys==15),refLabels==1),R((spdLabelsKeys==15),refLabels==1),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
p7=plot(L((spdLabelsKeys==16),refLabels==1),R((spdLabelsKeys==16),refLabels==1),'x','Color',daylightColours(7,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'Nature 5000K','Nature 5500K','Nature 6000K','Nature 6500K','Nature 7000K','Nature 7500K'});
xlabel('L/L+M');
ylabel('R/L+M');
%%

% there is probably a better way to do this but this will do for now
figure()
scatter3(L(:),S(:),I(:),'rx');
hold on;
scatter3(L(spdLabelsKeys>9),S(spdLabelsKeys>9),I(spdLabelsKeys>9),'bo');

% %%
% 
% for j=1:99
%     if refLabels(j)==1
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'kx');
%         hold on;
%         figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'kx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'kx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'kx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'kx');
%         hold on;
%     elseif refLabels(j)==2
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'bx');
%         hold on;
%         figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'bx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'bx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'bx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'bx');
%         hold on;
%     elseif refLabels(j)==3
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'rx');
%         hold on;
%                 figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'rx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'rx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'rx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'rx');
%         hold on;
%     elseif refLabels(j)==4
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'gx');
%         hold on;
%                 figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'gx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'gx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'gx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'gx');
%         hold on;
%     elseif refLabels(j)==5
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'cx');
%         hold on;
%                 figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'cx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'cx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'cx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'cx');
%         hold on;
%     elseif refLabels(j)==6
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'mx');
%         hold on;
%                 figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'mx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'mx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'mx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'mx');
%         hold on;
%     else
%         figure(3)
%         scatter3(L(:,j),S(:,j),I(:,j),'yx');
%         hold on;
%                 figure(4)
%         subplot(2,2,1)
%         plot(L(:,j),S(:,j),'yx');
%         hold on;
%         subplot(2,2,2)
%         plot(R(:,j),I(:,j),'yx');
%         hold on;
%         subplot(2,2,3)
%         plot(L(:,j),I(:,j),'yx');
%         hold on;
%         subplot(2,2,4)
%         plot(I(:,j),S(:,j),'yx');
%         hold on;
%     end
% end