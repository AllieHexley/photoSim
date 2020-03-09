% created by ACH 05/03/2020
% function to plot real world spectra in adapted MacLeod-Boynton
% (unlabelled)

% to do - add in covariance and fits to each plot

function [] = plotAllPhotoSim(L,S,I,R);

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

end