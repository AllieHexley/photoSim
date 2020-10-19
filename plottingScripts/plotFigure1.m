% plot Figure 1 of paper
% plots spectral sensitivities, example real-world specturm, example
% display primaries, example real-world signals, and example reproduced
% signals, with example introduced contrast
% created by ACH 01/07/2020

%% load data
ss = GetCIES026;
wlsCIES026 = (390:1:780)';
T_cies026 = ss(:,11:end);
T_cies026(isnan(T_cies026)) = 0;

%% set-up colors
lCol = [.5 0 0]; mCol = [0 0.5 0]; sCol = [0 0 .5];
rCol = [0 0.5 0.5]; iCol = [0.5 0 .5];

%% plot Fig1a - cone spectral sensitivities

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(390:780,T_cies026(1,:),'Color',sCol,'LineWidth',2);
h(2)=plot(390:780,T_cies026(2,:),'Color',mCol,'LineWidth',2);
h(3)=plot(390:780,T_cies026(3,:),'Color',lCol,'LineWidth',2);
legend(h,{'S','M','L'});
xlabel('Wavelength (nm)');
ylabel('Spectral Sensitivity');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1a.pdf','-dpdf');

%%  plot Fig1b - example real-world spectrum

fig = figure('defaultAxesFontSize',12);
% generate daylight spectrum
simRad = getSimulatedSpectra;
worldSpd = simRad(:,39230); % pick a random spectrum as an example from the simulated spectra
worldSpd = SplineSpd([390,5,79], worldSpd, [390,1,391]);
h(1)=plot(390:780,worldSpd,'k','LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1b.pdf','-dpdf');

%%  plot Fig1c - example three primary display primaries

fig = figure('defaultAxesFontSize',12);
% generate daylight spectrum
bbR = normpdf(390:780,450,(40./2.355));
bbG = normpdf(390:780,560,(40./2.355));
bbB = normpdf(390:780,630,(40./2.355));
wlsbb = [390:780];
spd = [bbR',bbG',bbB'];

% noramlise so area under primaries is 1
for i=1:size(spd,2)
    % calculate integral of illuminant spectra
    A(i) = trapz(wlsbb, spd(:,i));
    spd(:,i) = spd(:,i)./A(i);
end
hold on;
h(1)=plot(390:780,spd(:,1),'Color',sCol,'LineWidth',2);
h(2)=plot(390:780,spd(:,2),'Color',mCol,'LineWidth',2);
h(3)=plot(390:780,spd(:,3),'Color',lCol,'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1c.pdf','-dpdf');

%% calculate cone responses to worldSpd and display
lmsworldSpd = T_cies026(1:3,:)*worldSpd; % lms from worldSpd
LMS2RGB = T_cies026(1:3,:)*spd; % get RGB->LMS conversion for display
rgb = inv(LMS2RGB) * lmsworldSpd; % calculate RGB needed to match LMS to worldSpd
dispspd = rgb(1).*spd(:,1)+rgb(2).*spd(:,2)+rgb(3).*spd(:,3); % calculate output spectrum of display
lms3PDisp = T_cies026(1:3,:)*dispspd; % lms from display output
lmsri3PDisp = T_cies026*dispspd;
lmsriworldSpd = T_cies026*worldSpd;

%%  plot Fig1d - spectral output from display

fig = figure('defaultAxesFontSize',12);
h(1)=plot(390:780,dispspd,'Color',[0.6,0.6,0.6],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1d.pdf','-dpdf');

%% Fig 1e - cone responses to real-world spectrum and display spectrum

fig = figure('defaultAxesFontSize',12);
b = bar([lmsworldSpd'; lms3PDisp']','LineWidth',1.5);
b(1,1).FaceColor = 'flat';
b(1,1).EdgeColor = [0 0 0];
b(1,1).CData(3,:) = [.5 0 0];
b(1,1).CData(2,:) = [0 .5 0];
b(1,1).CData(1,:) = [0 0 .5];
b(1,2).FaceColor = 'flat';
b(1,2).EdgeColor = [0.5 0.5 0.5];
b(1,2).CData(3,:) = [.5 0 0];
b(1,2).CData(2,:) = [0 .5 0];
b(1,2).CData(1,:) = [0 0 .5];
xticklabels({'S','M','L'});
xlim([0.5,3.5]);
ylabel('Photoreceptor Response');
ylim([0,0.085]);
yticklabels({});
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1e.pdf','-dpdf');

%% Fig 1f - plot LMSRI spectral sensitivities

fig = figure('defaultAxesFontSize',12);
hold on;
h(1)=plot(390:780,T_cies026(4,:),'Color',rCol,'LineWidth',2);
h(2)=plot(390:780,T_cies026(5,:),'Color',iCol,'LineWidth',2);
legend(h,{'R','I'});
xlabel('Wavelength (nm)');
ylabel('Spectral Sensitivity');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1f.pdf','-dpdf');

%% calculate cone responses to real world spectrum and display
ri3PDisp = T_cies026(4:5,:)*dispspd;
riworldSpd = T_cies026(4:5,:)*worldSpd;

%% Fig 1g - rod and mel responses to real-world spectrum and display spectrum

fig = figure('defaultAxesFontSize',12);
b = bar([riworldSpd'; ri3PDisp']','LineWidth',1.5);
b(1,1).FaceColor = 'flat';
b(1,1).EdgeColor = [0 0 0];
b(1,1).CData(2,:) = iCol;
b(1,1).CData(1,:) = rCol;
b(1,2).FaceColor = 'flat';
b(1,2).EdgeColor = [0.6 0.6 0.6];
b(1,2).CData(2,:) = iCol;
b(1,2).CData(1,:) = rCol;
xticklabels({'R','I'});
xlim([0.5,2.5]);
ylabel('Photoreceptor Response');
ylim([0,0.085]);
yticklabels({});
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1g.pdf','-dpdf');

%%  plot Fig1h - example five primary display primaries

fig = figure('defaultAxesFontSize',12);
% generate daylight spectrum
bbR = normpdf(390:780,460,(40./2.355));
bbG = normpdf(390:780,530,(40./2.355));
bbB = normpdf(390:780,580,(40./2.355));
bbC = normpdf(390:780,610,(40./2.355));
bbM = normpdf(390:780,650,(40./2.355));
wlsbb = [390:780];
spd5 = [bbR',bbG',bbB',bbC',bbM'];

% noramlise so area under primaries is 1
for i=1:size(spd5,2)
    % calculate integral of illuminant spectra
    A(i) = trapz(wlsbb, spd5(:,i));
    spd5(:,i) = spd5(:,i)./A(i);
end
hold on;
h(1)=plot(390:780,spd5(:,1),'Color',sCol,'LineWidth',2);
h(2)=plot(390:780,spd5(:,2),'Color',mCol,'LineWidth',2);
h(3)=plot(390:780,spd5(:,5),'Color',lCol,'LineWidth',2);
h(3)=plot(390:780,spd5(:,4),'Color',[0.4940,0.1840,0.5560],'LineWidth',2);
h(3)=plot(390:780,spd5(:,3),'Color',[0.9290,0.6940,0.1250],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1h.pdf','-dpdf');

%% calculate cone responses to five primary display display
lmsriworldSpd = T_cies026*worldSpd; % lms from worldSpd
LMSRI2RGBCM = T_cies026*spd5; % get RGB->LMS conversion for display
rgbcm = inv(LMSRI2RGBCM) * lmsriworldSpd; % calculate RGB needed to match LMS to worldSpd
dispspd5 = rgbcm(1).*spd5(:,1)+rgbcm(2).*spd5(:,2)+rgbcm(3).*spd5(:,3)+rgbcm(4).*spd5(:,4)+rgbcm(5).*spd5(:,5); % calculate output spectrum of display
lmsri5PDisp = T_cies026*dispspd5; % lms from display output

%% plot Fig 1i - spectral output of 5P display

fig = figure('defaultAxesFontSize',12);
h(1)=plot(390:780,dispspd5,'Color',[0.25,0.25,0.25],'LineWidth',2);
xlabel('Wavelength (nm)');
ylabel('Relative Power');
yticklabels({});
xticks([400,500,600,700]);
xticklabels({'400','500','600','700'});
xlim([390,780]);
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1i.pdf','-dpdf');

%% plot Fig 1j - rod and mel responses to real-world spectrum and display spectrum

fig = figure('defaultAxesFontSize',12);
b = bar([lmsriworldSpd'; lmsri5PDisp']','LineWidth',1.5);
b(1,1).FaceColor = 'flat';
b(1,1).EdgeColor = [0 0 0];
b(1,1).CData(5,:) = iCol;
b(1,1).CData(4,:) = rCol;
b(1,1).CData(3,:) = lCol;
b(1,1).CData(2,:) = mCol;
b(1,1).CData(1,:) = sCol;
b(1,2).FaceColor = 'flat';
b(1,2).EdgeColor = [0.25 0.25 0.25];
b(1,2).CData(5,:) = iCol;
b(1,2).CData(4,:) = rCol;
b(1,2).CData(3,:) = lCol;
b(1,2).CData(2,:) = mCol;
b(1,2).CData(1,:) = sCol;
xticklabels({'S','M','L','R','I'});
xlim([0.5,5.5]);
ylabel('Photoreceptor Response');
ylim([0,0.085]);
yticklabels({});
axis square
grid on;
box on;
fig.PaperUnits = 'inches';
fig.PaperSize = [3.1,3.1];
fig.PaperPositionMode = 'manual';
fig.PaperPosition=[0.1 0.1 3 3];
print(fig, '..\plots\fig1j.pdf','-dpdf');