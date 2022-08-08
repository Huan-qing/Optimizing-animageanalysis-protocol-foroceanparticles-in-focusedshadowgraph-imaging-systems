

clear; close all; clc
% go to the path where the files is, and load the csv
% cd '/home/hhuan006/batch2021fall/T13/';
path='/home/hhuan006/datan2022batch/T30/sub10/';
list=dir(fullfile(path,'*cast-1_*'));
[~,index]=natsort({list.name});
list=list(index);
cat=cell(1,length(list));
m={};
% L = (zeros(0,34));
% L.Properties.VariableNames = {'image_number', 'k', 'Area', 'MajorAxisLength',  'MinorAxisLength', 'EquivDiameter', 'Eccentricity', 'Orientation', 'FilledArea', 'Solidity', 'Perimeter', 'MeanIntensity','date','time', 'Conductivity','Conductivity2','temperature', 'temperature2','pressure', 'oxy', 'oxy2','fluorescence[ug/l]',  'fluorescence[mg/m^3]', 'BeamAttenuation','BeamTransmission', 'PSU', 'PSU2','density','density2','depth','days','hours','min','sec'}
for t=1:length(list)
    
        name=list(t).name
        cat=readtable(strcat(path,name));
    if size(cat,1)~=0
        m=[m;cat];
    end
end
writetable(m,'/home/hhuan006/datan2022batch/T30/merge/date_cast-1_T30.csv');
