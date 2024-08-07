% function to find the photoreceptor signal distortions introduced when
% trying to reproduce the real-world spectra on a display
% created by ACH 01/07/2020

function [display] = getDistortions(matchedSignals,Sim,display,displayPrimaries,T_cies026,mb026,smallestBit,nPrimaries);

% inputs:
% 1) which photoreceotpr signals to match, vector including elements from 1-5 where 1=S;2=M;3=L;4=R;5=I
% 2) structure of simulated real world spectra
% 3) structure of display
% 4) primaries of the display
% 5) photoreceptor spectral sensitivities
% 6) MacLeod-Boynton chromaticity coordinates to produce extended
% MacLeod-Boynton diagram
% 7) smallest bit that display can achieve, for example, for an 8-bit display this would be 1./255
% 8) number of priamries in display

% outputs:
% display structure containing fields for distorted photoreceptor signals
% and which of those distorted signals can be reproduced on the display

% with constraint that we reproduce the signals specified in matchedSignals exactly, and let
% mel and rod distort freely

% calculate RGB to LMS conversion for the displays
if nPrimaries<=5
    RGB2SS = (T_cies026(matchedSignals,:)*displayPrimaries);
    SS2RGB = inv(RGB2SS);
elseif nPrimaries==6
    RGB2SS = (T_cies026(matchedSignals,:)*displayPrimaries(:,[1,2,3,5,6]));
    SS2RGB = inv(RGB2SS);
end

% for all real-world spectra LMS values, calculate mel and rod values
% you get when the real-world spectra are reproduced on the display
for i = 1:length(Sim.ss)
    if nPrimaries<=5
        % find out if the distorted photoreceptor signal can be reproduced
        % on the display with the conditons such that:
        % 1) No primary can be negative
        % 2) All of the primaries can't be 0 (as we allow intensity to
        % scale freely, it's not interesting if we can reproduce the real
        % world spectra at 0 intensity on the display)
        rgbSetting{i} = SS2RGB*Sim.ss(matchedSignals,i);
        rgbRealized{i} = (rgbSetting{i}'.*displayPrimaries); 
        if sum(rgbSetting{i}(:)<=0)>=1 | ...
                sum((rgbSetting{i}./max(rgbSetting{i}))<=(smallestBit./2))>=nPrimaries; % ratio required to produce one <1 gives the others at less than 1 bit
            display.ssReproducible(:,i) = 0;
        else
            display.ssReproducible(:,i)=1;
        end
        display.ssReproducible = logical(display.ssReproducible);
    else
        rgbSetting{i} = SS2RGB*Sim.ss(matchedSignals,i);
        rgbRealized{i} = (rgbSetting{i}'.*displayPrimaries(:,[1,2,3,5,6])); 
    end
    if nPrimaries == 3
        display.ssDistorted(:,i) = T_cies026*(rgbRealized{i}(:,1)+rgbRealized{i}(:,2)+rgbRealized{i}(:,3));
        % for MB too
        mb(:,i) = mb026*(rgbRealized{i}(:,1)+rgbRealized{i}(:,2)+rgbRealized{i}(:,3));
        display.mbDistorted(1,i) = mb(1,i)./(mb(2,i)+mb(3,i));
        display.mbDistorted(2,i) = mb(3,i)./(mb(2,i)+mb(3,i));
        display.mbDistorted(3,i) = mb(5,i)./(mb(2,i)+mb(3,i));
    elseif nPrimaries >=5
        display.ssDistorted(:,i) = T_cies026*(rgbRealized{i}(:,1)+rgbRealized{i}(:,2)+rgbRealized{i}(:,3)+rgbRealized{i}(:,4)+rgbRealized{i}(:,5));
        mb(:,i) = mb026*(rgbRealized{i}(:,1)+rgbRealized{i}(:,2)+rgbRealized{i}(:,3)+rgbRealized{i}(:,4)+rgbRealized{i}(:,5));
        display.mbDistorted(1,i) = mb(1,i)./(mb(2,i)+mb(3,i));
        display.mbDistorted(2,i) = mb(3,i)./(mb(2,i)+mb(3,i));
        display.mbDistorted(3,i) = mb(5,i)./(mb(2,i)+mb(3,i));
    end
end
   
    
    % for 6-primary displays iterate over all 6 primaries
    if nPrimaries==6
        C = nchoosek([1,2,3,4,5,6],5);
        reproducible = zeros(length(Sim.ss),6);
        for ii=1:length(C)
            spd = displayPrimaries(:,[C(ii,:)]);
            rgb2lms = T_cies026([1,2,3,4,5],:)*spd;
            lms2rgb = inv(rgb2lms);
            rgbcm = lms2rgb * Sim.ss([1,2,3,4,5],:);
            for i=1:length(Sim.ss)
                if sum(round(rgbcm(:,i),2)>=0)==5 | sum((rgbcm(:,i)./max(rgbcm(:,i))))<=smallestBit./2>=6;
                    reproducible(i,ii) = 1;
                end
            end
        end
        display.ssReproducible = sum(reproducible,2);
        for i=1:length(display.ssReproducible)
            if display.ssReproducible(i)>=1;
                display.ssReproducible(i)=1;
            end
        end
        display.ssReproducible = logical(display.ssReproducible);
        
    end
end

