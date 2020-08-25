% function to plot real-world chromaticities on CIE xy diagram along with CIE xy gamut of displays
% created by ACH 01/07/2020
      
%% load relevant data file

load('photosimPhotoreceptorDistortions_ReproduceLMS.mat');

%%  plot gamut of each display on CIE xy diagram
colsList = [0,0,0;0.75,0.75,0.75;0.5,0.5,0.5;0.25,0.25,0.25;0.9,0.9,0.9];
figure('defaultAxesFontSize',18)
plotChromaticity();
hold on;
h(1) = plotRealWorldChromaticity(Sim);
h(2)=plotChromaticityReproductionPerDisplay(CRT,1,colsList);
h(3)=plotChromaticityReproductionPerDisplay(LCD,2,colsList);
h(4)=plotChromaticityReproductionPerDisplay(DP,3,colsList);
%h(5)=plotChromaticityReproductionPerDisplay(FP1,4,colsList);
%h(6)=plotChromaticityReproductionPerDisplay(FP2,5,colsList);
xlabel('x');
ylabel('y');
l=legend(h,{'RealWorld','CRT','LCD','Display++'});
%l=legend(h,{'RealWorld','CRT','LCD','Display++','FP1','FP2'});
set(l,'color',[0.98,0.8,0.8]);
axis square

%% plot bar graph of % chromaticity reproducible on each display

figure('defaultAxesFontSize',18)
bar([CRT.chromaticityReproductionMetric;LCD.chromaticityReproductionMetric;DP.chromaticityReproductionMetric],'FaceColor',[0.5,0.5,0.5]);
xticks([1,2,3]);
xticklabels(['   CRT   ';'   LCD   ';'Display++']);
ylim([0,100]);
ylabel('Chromaticity Reproduction Metric (%)');
axis square

%% plotting functions

function h = plotChromaticityReproductionPerDisplay(display,d,colsList)

    h = plot(display.xyYMax(1,display.idx),display.xyYMax(2,display.idx),'Color',colsList(d,:),'LineWidth',2);

end

function g = plotRealWorldChromaticity(Sim)

    g = plot(Sim.xyY(1,:),Sim.xyY(2,:),'wo');

end
