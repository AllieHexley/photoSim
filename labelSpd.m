% created by ACH 05/03/2020
% function to label the illuminant spectra for the photosim project

function [spdLabels, spdLabelsKeys, spdLabelsValues] = labelSpd;

spdLabels = readtable('318Illuminants_sampleTypes.csv','ReadVariableNames',false);
spdLabels = table2cell(spdLabels);
% number labels
keySetSpd = {'Fluorescent Broadband', 'Fluorescent Narrowband', 'LED Phosphor', 'LED Mixed', 'LED Hybrid', 'High Intensity Discharge', 'Incandescent/Filament', 'Mathematical', 'Other'};
valueSetSpd = [1,2,3,4,5,6,7,8,9];
% define map
spdMap = containers.Map(keySetSpd,valueSetSpd);
%spdLabelsKeys = spdMap(spdLabels{1,1});

% add in additional labels for the Daylight spectra
% 10-D65 5000K; 11-D65 5500K, 12-D65 6000K, 13-D65 6500K, 14-D65 7000K, 15-D65 7500K 
for i=1:318
    spdLabelsKeys(i) = spdMap(spdLabels{1,i});
    spdLabelsValues{i} = spdLabels{1,i};
    if spdLabelsKeys(i) == 8;
        if strcmp('CIE D Series 5000 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 10;
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 5500 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 11;   
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 6000 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 12;   
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 6500 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 13;  
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 7000 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 14;
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 7500 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 15;
            spdLabelsValues{i} = spdLabels{2,i};
        elseif strcmp('CIE D Series 8000 K', spdLabels{2,i}) == 1
            spdLabelsKeys(i) = 16; 
            spdLabelsValues{i} = spdLabels{2,i};
        else
            continue;
        end
    end
end

end