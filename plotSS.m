% script to plot spectral sensitivities
% created by ACH 01/07/2020

%% load data
load('spectralSensitivities.mat');

%% 
figure('defaultAxesFontSize',18)
plotChromaticity();
xlabel('x');
ylabel('y');

%%  plot LMS
figure('defaultAxesFontSize',18)
hold on;
h(1)=plot(390:780,T_cies026(1,:),'b-','LineWidth',2);
h(2)=plot(390:780,T_cies026(2,:),'g-','LineWidth',2);
h(3)=plot(390:780,T_cies026(3,:),'r-','LineWidth',2);
legend(h,{'S','M','L'});
xlabel('Wavelength (nm)');
yticklabels({});
xlim([390,780]);
axis square

%%  plot LMSR
figure('defaultAxesFontSize',18)
hold on;
h(1)=plot(390:780,T_cies026(1,:),'b-','LineWidth',2,'Color',[0,0,1,0.1]);
h(2)=plot(390:780,T_cies026(2,:),'g-','LineWidth',2,'Color',[0,1,0,0.1]);
h(3)=plot(390:780,T_cies026(3,:),'r-','LineWidth',2,'Color',[1,0,0,0.1]);
h(4)=plot(390:780,T_cies026(4,:),'c-','LineWidth',2);
legend(h,{'S','M','L','R'});
xlabel('Wavelength (nm)');
yticklabels({});
xlim([390,780]);
axis square

%%  plot LMSRI
figure('defaultAxesFontSize',18)
hold on;
h(1)=plot(390:780,T_cies026(1,:),'b-','LineWidth',2,'Color',[0,0,1,0.1]);
h(2)=plot(390:780,T_cies026(2,:),'g-','LineWidth',2,'Color',[0,1,0,0.1]);
h(3)=plot(390:780,T_cies026(3,:),'r-','LineWidth',2,'Color',[1,0,0,0.1]);
h(4)=plot(390:780,T_cies026(4,:),'c-','LineWidth',2,'Color',[0,1,1,0.1]);
h(5)=plot(390:780,T_cies026(5,:),'m-','LineWidth',2);
legend(h,{'S','M','L','R','I'});
xlabel('Wavelength (nm)');
yticklabels({});
xlim([390,780]);
axis square