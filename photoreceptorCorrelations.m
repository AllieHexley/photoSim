% created by ACH 26/03/2020
% function to plot the correlations of L,M,S,R and I activity

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

%% reshape to single vector - rework this into function laters
% individual activations
L = reshape(photosim(:,:,3),[numObs,1]);
S = reshape(photosim(:,:,1),[numObs,1]);
I = reshape(photosim(:,:,5),[numObs,1]);
M = reshape(photosim(:,:,2),[numObs,1]);
R = reshape(photosim(:,:,4),[numObs,1]);
% photoreceptor responses matrix prm
prm = [S,M,L,R,I];

%% plot correlations

figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',14);
% histograms
subplot(5,5,13)
histogram(L,'FaceColor','r');
xlabel('L');
ylabel('');
subplot(5,5,7)
histogram(M,'FaceColor','g');
xlabel('M');
ylabel('');
subplot(5,5,1)
histogram(S,'FaceColor','b');
xlabel('S');
ylabel('');
subplot(5,5,19)
histogram(R,'FaceColor','c');
xlabel('R');
ylabel('');
subplot(5,5,25)
histogram(I,'FaceColor','m');
xlabel('I');
ylabel('');

% correlations for S
subplot(5,5,6)
plot(S,M,'b.');
xlabel('S');
ylabel('M');
[r,p]=corrcoef([S,M]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,2,txt);

subplot(5,5,11)
plot(S,L,'b.');
xlabel('S');
ylabel('L');
[r,p]=corrcoef([S,L]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,16)
plot(S,R,'b.');
xlabel('S');
ylabel('R');
[r,p]=corrcoef([S,R]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,21)
plot(S,I,'b.');
xlabel('S');
ylabel('I');
[r,p]=corrcoef([S,I]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

% for M

subplot(5,5,2)
plot(M,S,'g.');
xlabel('M');
ylabel('S');
[r,p]=corrcoef([S,M]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,12)
plot(M,L,'g.');
xlabel('M');
ylabel('L');
[r,p]=corrcoef([M,L]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,17)
plot(M,R,'g.');
xlabel('M');
ylabel('R');
[r,p]=corrcoef([M,R]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,22)
plot(M,I,'g.');
xlabel('M');
ylabel('I');
[r,p]=corrcoef([M,I]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

% for L

subplot(5,5,3)
plot(L,S,'r.');
xlabel('L');
ylabel('S');
[r,p]=corrcoef([L,S]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,8)
plot(L,M,'r.');
xlabel('L');
ylabel('M');
[r,p]=corrcoef([M,L]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,18)
plot(L,R,'r.');
xlabel('L');
ylabel('R');
[r,p]=corrcoef([L,R]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,23)
plot(L,I,'r.');
xlabel('L');
ylabel('I');
[r,p]=corrcoef([L,I]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

% for R

subplot(5,5,4)
plot(R,S,'c.');
xlabel('R');
ylabel('S');
[r,p]=corrcoef([R,S]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,9)
plot(R,M,'c.');
xlabel('R');
ylabel('M');
[r,p]=corrcoef([R,M]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,14)
plot(R,L,'c.');
xlabel('R');
ylabel('L');
[r,p]=corrcoef([L,R]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,24)
plot(R,I,'c.');
xlabel('R');
ylabel('I');
[r,p]=corrcoef([R,I]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

% for I

subplot(5,5,5)
plot(I,S,'m.');
xlabel('I');
ylabel('S');
[r,p]=corrcoef([I,S]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,10)
plot(I,M,'m.');
xlabel('I');
ylabel('M');
[r,p]=corrcoef([I,M]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,15)
plot(I,L,'m.');
xlabel('I');
ylabel('L');
[r,p]=corrcoef([I,L]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

subplot(5,5,20)
plot(I,R,'m.');
xlabel('I');
ylabel('R');
[r,p]=corrcoef([R,I]);
txt = {['R=' num2str(r(2))],['p=' num2str(p(2))]};
text(6,0.9,txt);

%% heatmap of correlations
[r,p]=corrcoef([S,M]);
rsm = r(2);
rss = r(1);
[r,p]=corrcoef([S,L]);
rsl = r(2);
[r,p]=corrcoef([S,R]);
rsr = r(2);
[r,p]=corrcoef([S,I]);
rsi = r(2);
[r,p]=corrcoef([M,L]);
rml = r(2);
rmm = r(1);
[r,p]=corrcoef([M,R]);
rmr = r(2);
[r,p]=corrcoef([M,I]);
rmi = r(2);
[r,p]=corrcoef([L,R]);
rlr = r(2);
rll = r(1);
[r,p]=corrcoef([L,I]);
rli = r(2);
[r,p]=corrcoef([R,I]);
rri = r(2);
rrr = r(1);
rii = r(1);

T = [rss,rsm,rsl,rsr,rsi;rsm,rmm,rml,rmr,rmi;rsl,rml,rll,rlr,rli;rsr,rmr,rlr,rrr,rri;rsi,rmi,rli,rri,rii];
% T.Properties.VariableNames = {'S','M','L','Rod','Mel'};
% T.Properties.RowNames = {'S','M','L','Rod','Mel'};
xvals = {'S','M','L','Rod','Mel'};
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',18);
h = heatmap(xvals,xvals,T,'FontSize',18);


