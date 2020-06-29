% percent in plots
% created by ACH 28/06/2020

%% plot xy diagrams for lum and mel and capture % of simulated in region

% define range of luminance levels for slice
lumLevs = [1,10,20,30,40,50,60,100,125,150,200,250];
melLevs = [0,0.02,0.04,0.06,0.08,0.1,0.12,0.15,0.2,0.25,0.3,0.325];

for i = i:length(lumLevs)
    [inLumCRT{i},onLumCRT{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),CRT.xyYFullRange(1,(CRT.xyYFullRange(3,:)>lumLevs(i))),CRT.xyYFullRange(2,(CRT.xyYFullRange(3,:)>lumLevs(i))));
    capturedLumCRT(i) = sum(inLumCRT{i})+sum(onLumCRT{i});
    percCapturedPerLumCRT(i) = round((capturedLumCRT(i)./totalSpec).*100,1);
    
    [inMelCRT{i},onMelCRT{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),CRT.xyYFullRange(1,(CRT.ssFullRange(5,:)>melLevs(i))),CRT.xyYFullRange(2,(CRT.ssFullRange(5,:)>melLevs(i))));
    capturedMelCRT(i) = sum(inMelCRT{i})+sum(onMelCRT{i});
    percCapturedPerMelCRT(i) = round((capturedMelCRT(i)./totalSpec).*100,1);
    
    [inLumLCD{i},onLumLCD{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),LCD.xyYFullRange(1,(LCD.xyYFullRange(3,:)>lumLevs(i))),LCD.xyYFullRange(2,(LCD.xyYFullRange(3,:)>lumLevs(i))));
    capturedLumLCD(i) = sum(inLumLCD{i})+sum(onLumLCD{i});
    percCapturedPerLumLCD(i) = round((capturedLumLCD(i)./totalSpec).*100,1);
    
    [inMelLCD{i},onMelLCD{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),LCD.xyYFullRange(1,(LCD.ssFullRange(5,:)>melLevs(i))),LCD.xyYFullRange(2,(LCD.ssFullRange(5,:)>melLevs(i))));
    capturedMelLCD(i) = sum(inMelLCD{i})+sum(onMelLCD{i});
    percCapturedPerMelLCD(i) = round((capturedMelLCD(i)./totalSpec).*100,1);
    
    [inLumDP{i},onLumDP{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),DP.xyYFullRange(1,(DP.xyYFullRange(3,:)>lumLevs(i))),DP.xyYFullRange(2,(DP.xyYFullRange(3,:)>lumLevs(i))));
    capturedLumDP(i) = sum(inLumDP{i})+sum(onLumDP{i});
    percCapturedPerLumDP(i) = round((capturedLumDP(i)./totalSpec).*100,1);
    
    [inMelDP{i},onMelDP{i}] = inpolygon(Sim.xyY(1,:),Sim.xyY(2,:),DP.xyYFullRange(1,(DP.ssFullRange(5,:)>melLevs(i))),DP.xyYFullRange(2,(DP.ssFullRange(5,:)>melLevs(i))));
    capturedMelDP(i) = sum(inMelDP{i})+sum(onMelDP{i});
    percCapturedPerMelDP(i) = round((capturedMelDP(i)./totalSpec).*100,1);
end

%%

save('inXY.mat','capturedLumCRT', 'capturedLumLCD', 'capturedLumDP', 'capturedMelCRT', 'capturedMelLCD', 'capturedMelDP', 'percCapturedPerLumCRT', 'percCapturedPerLumLCD', 'percCapturedPerLumDP', 'percCapturedPerMelCRT', 'percCapturedPerMelLCD', 'percCapturedPerMelDP');

%%
figure('defaultAxesFontSize',18)
sgtitle('Percentage Simulated Spectra Reproducible','FontSize',18);
subplot(1,2,1)
plot(lumLevs(2:end),percCapturedPerLumCRT(2:end),'k.-');
hold on;
plot(lumLevs(2:end),percCapturedPerLumLCD(2:end),'b.-');
plot(lumLevs(2:end),percCapturedPerLumDP(2:end),'r.-');
legend('CRT','LCD','DP')
xlabel('Lum');
ylabel('% Sim Captured');
set(gca,'XScale','log');
subplot(1,2,2)
plot(melLevs(2:end),percCapturedPerMelCRT(2:end),'k.-');
hold on;
plot(melLevs(2:end),percCapturedPerMelLCD(2:end),'b.-');
plot(melLevs(2:end),percCapturedPerMelDP(2:end),'r.-');
legend('CRT','LCD','DP')
xlabel('I');
ylabel('% Sim Captured');
set(gca,'XScale','log');