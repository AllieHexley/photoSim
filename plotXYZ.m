% plot XYZ
% created by ACH 01/07/2020

T_xyz = csvread('lin2012xyz10e_1_7sf.csv');
wls_xyz = T_xyz(:, 1);
T_xyz = 683*T_xyz(wls_xyz >= 390 & wls_xyz <= 780, 2:end)';
wls_xyz = wls_xyz(wls_xyz >= 390 & wls_xyz <= 780, 1);
wls_xyz = [];

figure('defaultAxesFontSize',18)
p(1)=plot(390:780,T_xyz(1,:),'r-','LineWidth',2);
hold on;
p(2)=plot(390:780,T_xyz(2,:),'g-','LineWidth',2);
p(3)=plot(390:780,T_xyz(3,:),'b-','LineWidth',2);
yticklabels({''});
legend(p,{'X','Y','Z'});
xlabel('Wavelength (nm)');
xlim([390,780]);

%% and plot RGB

[rgb,M] = XYZToSRGBPrimary(T_xyz);

figure('defaultAxesFontSize',18)
p(1)=plot(390:780,rgb(1,:),'r-','LineWidth',2);
hold on;
p(2)=plot(390:780,rgb(2,:),'g-','LineWidth',2);
p(3)=plot(390:780,rgb(3,:),'b-','LineWidth',2);
yticklabels({''});
legend(p,{'R','G','B'});
xlabel('Wavelength (nm)');
xlim([390,780]);