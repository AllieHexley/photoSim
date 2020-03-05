% created by ACH 05/03/2020
% function to label reflectance spectra for photosim project

function [refLabels, refLabelsKeys] = labelRef;

refLabels = csvread('99Reflectances_sampleTypes.csv');
refLabels = refLabels(1,:); %concatenate because reads in empty rows
%1-nature,2-skin,3-textiles,4-paints,5-plastic,6-printed,7-color system

% number labels
keySetRef = {'Nature','Skin','Textiles','Paints','Plastic','Printed','Color_System'};
valueSetRef = [1,2,3,4,5,6,7];
% define map
refMap=containers.Map(valueSetRef,keySetRef);

for j=1:99
    refLabelsKeys{j} = refMap(refLabels(j));
end

end