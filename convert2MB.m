% created by ACH 10/04/2020
% function to convert activations into MacLeod-Boynton Space
% normalized so that V(lambda)photopic; V(lambda) scotopic, and Smb and Imb peak at one

% load the photoreceptor spectral sensitivites
% [S,M,L,Rod,Mel]
[T_cies026, S_cies026] = GetCIES026;
wlsCIES026 = (390:1:780)';
% rescale only over range where we have cone fundmanetals i.e. 390nm:780nm
% scale for spds with 1nm spacing
T_cies026_1nm = T_cies026(:,11:end);
% remove Nans
T_cies026_1nm(isnan(T_cies026_1nm)) = 0;

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
figure
plot(wlsCIES026, vLambda);
hold on;
plot(wlsCIES026,vLamSS(1:391,2),'b--');
title('VLambda check');

% now calculate chromaticity coordinates
lmb = (cl.*l)./vLambda;
mmb = (cm.*m)./vLambda;
smb = (cs.*s)./vLambda;

% plot chromaticity coordinates
figure
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

% plot chromaticity coordinates
figure
plot(lCoord, sCoord, 'k');
title('Traditional MacLeod Boynton Space');

