% created by ACH 05/03/2020
% function to label the illuminant spectra for the photosim project

function [spdLabels, spdLabelsValues] = labelSpd;

spdLabels = readtable('401Illuminants_sampleTypes.csv','ReadVariableNames',false);
spdLabels = table2cell(spdLabels);
% number labels
keySetSpd = {'LP-R','LM-R','FB-R','FN-R','HI-R','TF-R','LP-T','LM-T','FL-T','BB-T','DS-T','OT-T'};
valueSetSpd = {'LED Phosphor Real', 'LED Mixed Real', 'Fluorescent Broadband',' Fluorescent Narrowband', 'HID', 'Tungsten Filament', 'LED Phosphor Models' ,'LED Mixed Models', 'Fluorescent Models', 'Blabkbody Radiation' ,'D-Series Illuminant', 'Other'};
% define map
spdMap = containers.Map(keySetSpd,valueSetSpd);
%spdLabelsKeys = spdMap(spdLabels{1,1});

% add in additional labels for the Daylight spectra
% 10-D65 5000K; 11-D65 5500K, 12-D65 6000K, 13-D65 6500K, 14-D65 7000K, 15-D65 7500K 
for i=1:401
    spdLabelsValues{i} = spdMap(spdLabels{i});
%     if spdLabels{i} == 'DS-T';
%         if strcmp('CIE D Series 5000 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 10;
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 5500 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 11;   
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 6000 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 12;   
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 6500 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 13;  
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 7000 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 14;
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 7500 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 15;
%             spdLabelsValues{i} = spdLabels{2,i};
%         elseif strcmp('CIE D Series 8000 K', spdLabels{2,i}) == 1
%             spdLabelsKeys(i) = 16; 
%             spdLabelsValues{i} = spdLabels{2,i};
%         else
%             continue;
%         end
%     end
end

end