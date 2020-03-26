% created by ACH 26/03/2020
% function to plot the spectral locus in adapted MacLeod Boynton space

%% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (390:1:780)';
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
T_cies026 = T_cies026(:,11:end);
% remove Nans
T_cies026(isnan(T_cies026)) = 0;

%% simulate spectral locus
slSPD = zeros(391,391);
for ii=1:391
    slSPD(ii,ii) = 1;
end;

% calculate spectral locus in L,M,S,R,I coordinate space
for m=1:391
    for k=1:5
        spectralLocus(m,k)=T_cies026(k,:)*(slSPD(:,m));
    end
end
% remove Nans
spectralLocus(isnan(spectralLocus))=0;

%% calculate macleod boynton of spectral locus
% for all need to shift so starts at 390 or else will reach infintiy
for n=1:391
    Lsl(n) = spectralLocus(n,3)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Msl(n) = spectralLocus(n,2)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Ssl(n) = spectralLocus(n,1)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Rsl(n) = spectralLocus(n,4)./((spectralLocus(n,3)+spectralLocus(n,2)));
    Isl(n) = spectralLocus(n,5)./((spectralLocus(n,3)+spectralLocus(n,2)));
end

%% normalise S, R, and I to unity
Ssl = Ssl./(max(Ssl));
Rsl = Rsl./(max(Rsl));
Isl = Isl./(max(Isl));

%% plot

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