% created by ACH 26/03/2020
% function to plot the spectral and daylight locus in adapted MacLeod
% Boynton space

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

%% plot spectral and daylight locus for just coorrelation space i.e. not MacLeod Boynton

% colours for daylight plot
daylightColours = [0.2,1,1;0.2,0.8,1;0.2,0.6,1.0;0.2,0.4,1.0;0.2,0.2,1.0;0.2,0.0,1.0];

%% plot spectral locus alone
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',12);
subplot(2,2,1)
plot(Lsl(:),Ssl(:),'r-');
hold on;
plot(Lsl(1:10:end),Ssl(1:10:end),'k.');
% add in some spectral loci
%wavs = {'390nm','400nm','410nm','420nm','430nm','440nm','450nm','460nm','470nm','480nm','490nm','500nm','510nm','520nm','530nm','540nm','550nm','560nm','570nm','580nm','590nm','600nm','610nm','620nm','630nm','640nm','650nm','660nm','670nm','680nm','690nm','700nm','710nm','720nm','730nm','740nm','750nm','760nm','770nm','780nm'};
wavs = {'400nm','450nm','500nm','550nm','600nm','650nm'};
c=1;
for i=11:50:261;
    text(Lsl(i),Ssl(i),wavs{c});
    c=c+1;
end
hold on;
p1=plot(Ldl(1),Sdl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Ldl(2),Sdl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Ldl(3),Sdl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Ldl(4),Sdl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Ldl(5),Sdl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Ldl(6),Sdl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,2)
plot(Rsl(:),Isl(:),'r-');
hold on;
plot(Rsl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Rsl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot(Rdl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Rdl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Rdl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Rdl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Rdl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Rdl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('R/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,3)
plot(Lsl(:),Isl(:),'r-');
hold on;
plot(Lsl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Lsl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot(Ldl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Ldl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Ldl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Ldl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Ldl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Ldl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L/L+M');
ylabel('Mel/L+M');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,4)
plot(Isl(:),Ssl(:),'r-');
hold on;
plot(Isl(1:10:end),Ssl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Isl(i),Ssl(i),wavs{c});
    c=c+1;
end
p1=plot(Idl(1),Sdl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Idl(2),Sdl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Idl(3),Sdl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Idl(4),Sdl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Idl(5),Sdl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Idl(6),Sdl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('Mel/L+M');
ylabel('S/L+M');

%% plot adapted MacLeod Boynton in 3D
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',12);
plot3(Lsl(:),Ssl(:),Isl(:),'r-');
hold on;
plot3(Lsl(1:10:end),Ssl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Lsl(i),Ssl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot3(Ldl(1),Sdl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot3(Ldl(2),Sdl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot3(Ldl(3),Sdl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot3(Ldl(4),Sdl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot3(Ldl(5),Sdl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot3(Ldl(6),Sdl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L/L+M');
ylabel('S/L+M');
zlabel('Mel/L+M');

%% and for just correlations
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
    Lsl(n) = spectralLocus(n,3);
    Msl(n) = spectralLocus(n,2);
    Ssl(n) = spectralLocus(n,1);
    Rsl(n) = spectralLocus(n,4);
    Isl(n) = spectralLocus(n,5);
end

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
    Ldl(n) = daylightLocus(n,3);
    Mdl(n) = daylightLocus(n,2);
    Sdl(n) = daylightLocus(n,1);
    Rdl(n) = daylightLocus(n,4);
    Idl(n) = daylightLocus(n,5);
end

Ldl = Ldl./max(Ldl);
Mdl = Mdl./max(Mdl);
Sdl = Sdl./max(Sdl);
Rdl = Rdl./max(Rdl);
Idl = Idl./max(Idl);

%% plot spectral and daylight locus for just coorrelation space i.e. not MacLeod Boynton

% colours for daylight plot
daylightColours = [0.2,1,1;0.2,0.8,1;0.2,0.6,1.0;0.2,0.4,1.0;0.2,0.2,1.0;0.2,0.0,1.0];

%% plot spectral locus alone
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',12);
subplot(2,2,1)
plot(Lsl(:),Ssl(:),'r-');
hold on;
plot(Lsl(1:10:end),Ssl(1:10:end),'k.');
% add in some spectral loci
%wavs = {'390nm','400nm','410nm','420nm','430nm','440nm','450nm','460nm','470nm','480nm','490nm','500nm','510nm','520nm','530nm','540nm','550nm','560nm','570nm','580nm','590nm','600nm','610nm','620nm','630nm','640nm','650nm','660nm','670nm','680nm','690nm','700nm','710nm','720nm','730nm','740nm','750nm','760nm','770nm','780nm'};
wavs = {'400nm','450nm','500nm','550nm','600nm','650nm','700nm','750nm'};
c=1;
for i=11:50:261;
    text(Lsl(i),Ssl(i),wavs{c});
    c=c+1;
end
hold on;
p1=plot(Ldl(1),Sdl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Ldl(2),Sdl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Ldl(3),Sdl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Ldl(4),Sdl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Ldl(5),Sdl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Ldl(6),Sdl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L');
ylabel('S');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,2)
plot(Rsl(:),Isl(:),'r-');
hold on;
plot(Rsl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Rsl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot(Rdl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Rdl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Rdl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Rdl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Rdl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Rdl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('R');
ylabel('Mel');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,3)
plot(Lsl(:),Isl(:),'r-');
hold on;
plot(Lsl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Lsl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot(Ldl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Ldl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Ldl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Ldl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Ldl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Ldl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L');
ylabel('Mel');

%% plot adapted MacLeod Boynton space in rod and mel space
subplot(2,2,4)
plot(Isl(:),Ssl(:),'r-');
hold on;
plot(Isl(1:10:end),Ssl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Isl(i),Ssl(i),wavs{c});
    c=c+1;
end
p1=plot(Idl(1),Sdl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot(Idl(2),Sdl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot(Idl(3),Sdl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot(Idl(4),Sdl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot(Idl(5),Sdl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot(Idl(6),Sdl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('Mel');
ylabel('S');

%% plot adapted MacLeod Boynton in 3D
figure('DefaultAxesFontSize',18)
set(0,'DefaultTextFontSize',12);
plot3(Lsl(:),Ssl(:),Isl(:),'r-');
hold on;
plot3(Lsl(1:10:end),Ssl(1:10:end),Isl(1:10:end),'k.');
c=1;
for i=11:50:261;
    text(Lsl(i),Ssl(i),Isl(i),wavs{c});
    c=c+1;
end
p1=plot3(Ldl(1),Sdl(1),Idl(1),'x','Color',daylightColours(1,:),'MarkerSize',8,'LineWidth',2);
p2=plot3(Ldl(2),Sdl(2),Idl(2),'x','Color',daylightColours(2,:),'MarkerSize',8,'LineWidth',2);
p3=plot3(Ldl(3),Sdl(3),Idl(3),'x','Color',daylightColours(3,:),'MarkerSize',8,'LineWidth',2);
p4=plot3(Ldl(4),Sdl(4),Idl(4),'x','Color',daylightColours(4,:),'MarkerSize',8,'LineWidth',2);
p5=plot3(Ldl(5),Sdl(5),Idl(5),'x','Color',daylightColours(5,:),'MarkerSize',8,'LineWidth',2);
p6=plot3(Ldl(6),Sdl(6),Idl(6),'x','Color',daylightColours(6,:),'MarkerSize',8,'LineWidth',2);
legend([p1,p2,p3,p4,p5,p6],{'5000K','5500K','6000K','7000K','7500K','8000K'});
xlabel('L');
ylabel('S');
zlabel('Mel');
