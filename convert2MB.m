% created by ACH 10/04/2020
% function to convert activations into MacLeod-Boynton Space

% throughout rename activations to sensitivities

% normalized so that V(lambda)photopic; V(lambda) scotopic, and Smb and Imb peak at one
clear all;
clc;
% load the photoreceptor spectral sensitivites
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
activations = T_cies026_5nm;

% write out l,m,s,r,i for convenience
l=T_cies026_1nm(3,:);
m=T_cies026_1nm(2,:);
s=T_cies026_1nm(1,:);
r=T_cies026_1nm(4,:);
i=T_cies026_1nm(5,:);

% constants
cl = 0.692839;
cm = 0.349676;
cs = 0.0554786;

% calculate VLambda
vLambda = (cl.*l)+(cm.*m);

% check its the same as stockmon sharpe VLambda
vLamSS = csvread('linCIE2008v10e_1.csv');

% it is!

% plot VLambda
% figure
% plot(wlsCIES026, vLambda);
% hold on;
% plot(wlsCIES026,vLamSS(1:391,2),'b--');
% title('VLambda check');

% now calculate chromaticity coordinates
lmb = (cl.*l)./vLambda;
mmb = (cm.*m)./vLambda;
%smb = (cs.*s)./vLambda;
smb = s./vLambda;
smb = smb./max(smb);

% plot chromaticity coordinates
figure('defaultAxesFontSize',18)
plot(wlsCIES026, lmb, 'r');
hold on;
plot(wlsCIES026, mmb, 'g');
plot(wlsCIES026, smb, 'b');
title('chromaticity coordinates check');
mbSS = csvread('mb10_1.csv');
plot(wlsCIES026, mbSS(1:391,2), 'k--');
hold on;
plot(wlsCIES026, mbSS(1:391,3), 'k--');
plot(wlsCIES026, mbSS(1:391,4), 'k--');

% they match!

% now calculate MacLeod Boynton space i.e. L/L+M and S/L+M
lCoord = lmb./(lmb+mmb);
sCoord = smb./(lmb+mmb);

% try to caluclate melanopsin coordinate
imb = i./vLambda;
% calculate melanopsin constant
ci = 1./max(imb);
imb = imb./(max(imb));
plot(wlsCIES026, imb, 'm');

% calcualte melanopsin coord
iCoord = imb./(lmb+mmb);

% plot i vs V(lambda)
figure('defaultAxesFontSize',18)
plot(vLambda,i,'r-');
xlabel('L+M');
title('Spectral locus Mel vs Luminance');
ylabel('Mel');

%% plot chromaticity coordinates
figure('defaultAxesFontSize',18)
subplot(1,2,1);
plot(lCoord, sCoord, 'k');
title('Traditional MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('S/L+M');
subplot(1,2,2)
plot(lCoord, iCoord, 'k');
title('Adapted MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('I/L+M');

figure('defaultAxesFontSize',18)
plot3(lCoord,sCoord,iCoord,'k-')
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');

% get vector of mb coords
mbCoords = [smb(1,1:5:391); mmb(1,1:5:391); lmb(1,1:5:391); imb(1,1:5:391)];

%% get daylight locus
% simulate daylight locus
daylightSPD = csvread('daylightSpectra.csv');
wlsSpd = daylightSPD(3:end,1);
% cols are different illuminants, rows are diff wavelengths
daylightSPD = daylightSPD(3:end,2:end);

% normalise daylight locus
normSpd = normIllSpd(daylightSPD, wlsSpd);

% scale sensitivities appropriately
sensitivities = [activations(1,:).*cs;activations(2,:).*cm;activations(3,:).*cl;activations(5,:).*ci];

% this method works - do this for simulations as well
% calculate daylight locus in L,M,S,R,I coordinate space
for m=1:size(daylightSPD,2)
    for k=1:size(sensitivities,1)
        daylightLocus(m,k)=sensitivities(k,:)*(normSpd(:,m));
        daylightLocusMB(m,k) = mbCoords(k,:)*normSpd(:,m);
    end
end
% remove Nans
daylightLocus(isnan(daylightLocus))=0;

% get MacLeod Boynton of daylight Locus
L = daylightLocus(:,3)./ (daylightLocus(:,3)+daylightLocus(:,2));
S = daylightLocus(:,1)./ (daylightLocus(:,3)+daylightLocus(:,2));
I = daylightLocus(:,4)./ (daylightLocus(:,3)+daylightLocus(:,2));
% plot MacLeod Boynton of daylight locus
%% plot chromaticity coordinates
% REDO don't want to multiple by MC coord but want to convert activation to
% MB space

figure('defaultAxesFontSize',18)
subplot(1,2,1);
plot(lCoord, sCoord, 'r');
hold on;
plot(L,S,'bx');
title('Traditional MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('S/L+M');
legend('Spectral Locus', 'Daylight Locus');
subplot(1,2,2)
plot(lCoord, iCoord, 'r');
hold on;
plot(L,I,'bx');
title('Adapted MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('I/L+M');
legend('Spectral Locus', 'Daylight Locus');

figure('defaultAxesFontSize',18)
plot3(lCoord,sCoord,iCoord,'r-');
hold on;
plot3(L,S,I,'bx');
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');
legend('Spectral Locus', 'Daylight Locus');
%% get it for all simulations

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

% calculate L,M,S,R,I responses
for i=1:size(normSpd,2)
    for j=1:size(ref,2)
        for k=1:size(sensitivities,1)
            photosim(i,j,k)=sensitivities(k,:)*(ref(:,j).*normSpd(:,i));
        end
    end
end
% remove Nans
photosim(isnan(photosim))=0;

% get MacLeod Boynton of daylight Locus
Lp = photosim(:,:,3)./ (photosim(:,:,3)+photosim(:,:,2));
Sp = photosim(:,:,1)./ (photosim(:,:,3)+photosim(:,:,2));
Ip = photosim(:,:,4)./ (photosim(:,:,3)+photosim(:,:,2));
% plot MacLeod Boynton of daylight locus
%% plot chromaticity coordinates
% REDO don't want to multiple by MC coord but want to convert activation to
% MB space

figure('defaultAxesFontSize',18)
subplot(1,2,1);
hold on;
p1=plot(Lp,Sp,'k.');
p2=plot(L,S,'bx');
p3=plot(lCoord, sCoord, 'r');
title('Traditional MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('S/L+M');
legend([p1;p2;p3],{'Simulated Spectra'; 'Daylight Locus   ';'Spectral Locus   '});
subplot(1,2,2)
hold on;
p1=plot(Lp,Ip,'k.');
p2=plot(L,I,'bx');
p3=plot(lCoord, iCoord, 'r');
title('Adapted MacLeod Boynton Space');
xlabel('L/L+M');
ylabel('I/L+M');
legend([p1;p2;p3],{'Simulated Spectra'; 'Daylight Locus';'Spectral Locus'});

figure('defaultAxesFontSize',18)
scatter3(reshape(Lp,[39699,1]),reshape(Sp,[39699,1]),reshape(Ip,[39699,1]),'k.');
hold on;
plot3(lCoord,sCoord,iCoord,'r-')
plot3(L,S,I,'bx');
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('I/L+M');
legend('Simulated Spectra', 'Daylight Locus','Spectral Locus');
% work out the CIE melanopic efficacy of luminance stuff and check getting
% ball park values

% do hyperspectral images first though

% %% now try to calculate rod coord assume rod coord is scotopic V(lambda) and check this later
% figure()
% plot(wlsCIES026,r)
% 
% % luminous effiiency constants in lm/W
% photopic = 683;
% scotopic = 1700;
% 
% % scale so scotopic function is ratio greater than photopic
% rmb = (scotopic./photopic);
% 
% figure(5)
% plot(wlsCIES026, rmb, 'k')