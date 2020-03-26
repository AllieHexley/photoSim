% created by ACH 26/03/2020
% function to plot the spectral locus in adapted MacLeod Boynton space

%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (390:5:780)';
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
T_cies026 = T_cies026(:,11:5:end);
% remove Nans
T_cies026(isnan(T_cies026)) = 0;

%% simulate daylight locus
daylightSPD = csvread('daylightSpectra.csv');
wlsSpd = daylightSPD(3:end,1);
% cols are different illuminants, rows are diff wavelengths
daylightSPD = daylightSPD(3:end,2:end);

%% calculate daylight locus in L,M,S,R,I coordinate space
for m=1:6
    for k=1:5
        daylightLocus(m,k)=T_cies026(k,:)*(daylightSPD(:,m));
    end
end
% remove Nans
daylightLocus(isnan(daylightLocus))=0;

%% calculate macleod boynton of spectral locus
% for all need to shift so starts at 390 or else will reach infintiy
for n=1:6
    Ldl(n) = daylightLocus(n,3)./((daylightLocus(n,3)+daylightLocus(n,2)));
    Mdl(n) = daylightLocus(n,2)./((daylightLocus(n,3)+daylightLocus(n,2)));
    Sdl(n) = daylightLocus(n,1)./((daylightLocus(n,3)+daylightLocus(n,2)));
    Rdl(n) = daylightLocus(n,4)./((daylightLocus(n,3)+daylightLocus(n,2)));
    Idl(n) = daylightLocus(n,5)./((daylightLocus(n,3)+daylightLocus(n,2)));
end

%% plot daylight alone
figure()
subplot(2,2,1)
plot(Ldl(:),Sdl(:),'ro');
xlabel('L/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,2)
plot(Rdl(:),Idl(:),'ro');
xlabel('R/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,3)
plot(Ldl(:),Idl(:),'ro');
xlabel('L/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,4)
plot(Idl(:),Sdl(:),'ro');
xlabel('Mel/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton in 3D
figure()
scatter3(Ldl(:),Sdl(:),Idl(:),'ro');
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('Mel/L+M');