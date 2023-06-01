

 



tic

% date_time data
    clear all;
    close all;
    myFolder = ('./RedSea_analysis/diff/cast-5/');
    if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
    end
    % Change to whatever pattern you need.
    filePattern = fullfile(myFolder, '*.tif*'); 
    theFiles = dir(filePattern);
    [~, Index] = natsort({theFiles.name});
    theFiles   = theFiles(Index);
    n=0;
    M={};
    sec_total=500;
    sec=454;
    each=round(length(theFiles)/sec_total);
    savepath='./RedSea_analysis_test/batch/SPACE/T=29/sub500/cast-5/';
    if ~exist(savepath)
      mkdir(savepath)
    end
     % to 11
    if (sec-1)*each>=length(theFiles)
        writetable(cell2table(M),[savepath 'T=29_date_time-cast-5_sec=454.csv']); 
    else
        if sec*each<length(theFiles)
            last_file=each*sec;
        else
        last_file=length(theFiles);
        end
    
for i=each*(sec-1)+1:last_file
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s', fullFileName);
    

    % set the threshold for binary
    T=29;
    % get the path, name(mainly we get) and file extensions
    [filepath,name,ext] = fileparts(baseFileName);
    % seperate the data based on the ' ' into image number, date and time
    time_ctd = regexp(name,' ','split');
    % get the variable of image number
    image_label = time_ctd{1};
    % read the number after image
    [image_number,~,~]=sscanf(image_label,'diff_gain%f');
    % read the date on the name of image
    % read the date on the name of image
    dat_{i}=baseFileName(end-21:end-13);
    % read the time on the name of image
    tim_{i}=baseFileName(end-12:end-5);
    % change the '-' to ':' on name of time
    tim_{i}= strrep(tim_{i},'-',':');
    new=baseFileName(end-21:end);
    folder = fileparts(which(baseFileName)); % Determine where demo folder is (works with all versions).
    fullFileName = fullfile(folder, baseFileName);
    originalImage = imread(strcat(myFolder,fullFileName));
    % Check to make sure that it is grayscale, just in case the user substituted their own image.
    [rows, columns, numberOfColorChannels] = size(originalImage);
    thresholdValue = T;
    binaryImage = rem_simcan(originalImage,  thresholdValue); 
%    binaryImage = originalImage>thresholdValue; 
%    imshow(binaryImage)

blobMeasurements = regionprops(binaryImage, originalImage, 'all');%ABB
numberOfBlobs = size(blobMeasurements, 1);

%fprintf(1,'Image% Blob%      Mean Intensity  Area   Perimeter    Centroid       Diameter');
%M=cell(numberOfBlobs,14); %preload matrix M with zeros
% Loop over all blobs printing their measurements to the command window.
for k = 1 : numberOfBlobs           % Loop through all blobs.
%   fprintf(1,'%2d %2d %17.1f %11.1f %8.1f %8.1f %8.1f % 8.1f', image_number, k, meanGL, blobArea, blobPerimeter, blobCentroid, blobECD(k));%ABB
% m=zeros(row_num,column_num);
% m(i,:) = new_row;
% M=[];
Area = blobMeasurements(k).Area;
MajorAxisLength=blobMeasurements(k).MajorAxisLength;
MinorAxisLength=blobMeasurements(k).MinorAxisLength;
EquivDiameter=blobMeasurements(k).EquivDiameter;
Eccentricity=blobMeasurements(k).Eccentricity;
Orientation=blobMeasurements(k).Orientation;
FilledArea=blobMeasurements(k).FilledArea;
Solidity=blobMeasurements(k).Solidity;
Perimeter=blobMeasurements(k).Perimeter;
MeanIntensity=blobMeasurements(k).MeanIntensity;
Date = dat_{i};
Timepoint = tim_{i};
% save all the measurement results to a new cell
[M{k+n,1}]=[image_number];
[M{k+n,2}]=[k];
[M{k+n,3}]=[Area];
[M{k+n,4}]=[MajorAxisLength];
[M{k+n,5}]=[MinorAxisLength];
[M{k+n,6}]=[EquivDiameter ];
[M{k+n,7}]=[Eccentricity];
[M{k+n,8}]=[ Orientation];
[M{k+n,9}]=[ FilledArea ];
[M{k+n,10}]=[Solidity];
[M{k+n,11}]=[Perimeter];
[M{k+n,12}]=[MeanIntensity];
[M{k+n,13}]=[Date];
[M{k+n,14}]=[Timepoint];
end
% if there are no particles in the image, the values would be 0
if isempty(k)==1
Date = dat_{i};
Timepoint = tim_{i};
k=1;
[M{k+n,1}]=[image_number];
[M{k+n,2}]=[0];
[M{k+n,3}]=[0];
[M{k+n,4}]=[0];
[M{k+n,5}]=[0];
[M{k+n,6}]=[0];
[M{k+n,7}]=[0];
[M{k+n,8}]=[0];
[M{k+n,9}]=[0];
[M{k+n,10}]=[0];
[M{k+n,11}]=[0];
[M{k+n,12}]=[0];
[M{k+n,13}]=[Date];
[M{k+n,14}]=[Timepoint];
end
n=n+k;
end
M=sortrows(M,1);
L=cell2table(M);

writetable(L,[savepath 'T=29_date_time-cast-5_sec=454.csv']); 
end




%%





clear; close all; clc
% go to the path where the files is, and load the csv
path= './RedSea_analysis_test/batch/SPACE/T=29/sub500/cast-5/';
path_ctd='/home/hhuan006/CTD_redsea/';
ctd_list=dir(fullfile(path_ctd,'*cast-5.*'));
ctd_name=ctd_list(1).name;
data_ctd=readcell(strcat(path_ctd,ctd_name));
data_ctd=data_ctd(2:end,:);
filelist=dir(fullfile(path,'*cast-5_sec=454.csv'));
n=1;
Comb={};
m=1;
L={};
% L = (zeros(0,34));
% L.Properties.VariableNames = {'image_number', 'k', 'Area', 'MajorAxisLength',  'MinorAxisLength', 'EquivDiameter', 'Eccentricity', 'Orientation', 'FilledArea', 'Solidity', 'Perimeter', 'MeanIntensity','date','time', 'time_point','pressure','temperature', 'salinity','fluorescent', 'oxy', 'PAR','potential temp',  'density', 'depth','OBS', 'beam attenuation', 'transmission','voltage','hours','min','sec','julian days'}
ctd_time=string(data_ctd(:,1));
for i = 1: length(ctd_time)
    timepoint_ctd = regexp(ctd_time(i),':','split');
    if ~ismissing(timepoint_ctd)
    hour_ctd(i,1)=str2num(timepoint_ctd(1));
    minute_ctd(i,1)=str2num(timepoint_ctd(2));
    second_ctd(i,1)=str2num(timepoint_ctd(3));
    else
    hour_ctd(i,1)=NaN;
    minute_ctd(i,1)=NaN;
    second_ctd(i,1)=NaN;
    end
end
file_name=filelist(1).name;    
data_particle=readtable(strcat(path,file_name));
save_name=strrep(file_name,'date_time','combine');
% get the time in image data
time_point=table2cell(data_particle(:,14));
[w,p]=size(time_point);
% change the format of x into string as B
time_p=string(time_point);
for j=1:length(time_point)

% get the row number and column number of the CTD data


% give the data of Julian days into a variable A

% calculate the julian days into hour, minute and second, and save into
% seperate variables

% split the string in the image data into hour, minute and second
unt=time_p(j,1);
P=strsplit(unt,':');
hour=str2num(P{1});
hour=hour;

minute=str2num(P{2});
second=str2num(P{3});
% if the time points in 2 variabls match, the write down the data
if ~isempty(find(hour_ctd==hour&minute_ctd==minute&second_ctd==second))==1
        i=find(hour_ctd==hour&minute_ctd==minute&second_ctd==second); i=i(1);
 
        Comb{n,1}=table2cell(data_particle(j,1));
        Comb{n,2}=table2cell(data_particle(j,2));
        Comb{n,3}=table2cell(data_particle(j,3));
        Comb{n,4}=table2cell(data_particle(j,4));
        Comb{n,5}=table2cell(data_particle(j,5));
        Comb{n,6}=table2cell(data_particle(j,6));
        Comb{n,7}=table2cell(data_particle(j,7));
        Comb{n,8}=table2cell(data_particle(j,8));
        Comb{n,9}=table2cell(data_particle(j,9));
        Comb{n,10}=table2cell(data_particle(j,10));
        Comb{n,11}=table2cell(data_particle(j,11));
        Comb{n,12}=table2cell(data_particle(j,12));
        Comb{n,13}=table2cell(data_particle(j,13));
        Comb{n,14}=table2cell(data_particle(j,14));
        Comb{n,15}=(data_ctd(i,1));
        Comb{n,16}=(data_ctd(i,2));
        Comb{n,17}=(data_ctd(i,3));
        Comb{n,18}=(data_ctd(i,4));
        Comb{n,19}=(data_ctd(i,5));
        Comb{n,20}=(data_ctd(i,6));
        Comb{n,21}=(data_ctd(i,7));
        Comb{n,22}=(data_ctd(i,8));
        Comb{n,23}=(data_ctd(i,9));
        Comb{n,24}=(data_ctd(i,10));
        Comb{n,25}=(data_ctd(i,11));
        Comb{n,26}=(data_ctd(i,12));
        Comb{n,27}=(data_ctd(i,13));
        Comb{n,28}=(data_ctd(i,14));
        Comb{n,29}=(data_ctd(i,15));

    
  
        %[image_number k Area MajorAxisLength  MinorAxisLength EquivDiameter Eccentricity Orientation FilledArea Solidity Perimeter MeanIntensity data time Conductivity   Conductivity2  temperature temperature2   pressure    oxy oxy2   fluorescence[ug/l]  fluorescence[mg/m^3]  BeamAttenuation  BeamTransmission PSU PSU2   density density2   depth   days    hours   min  sec ]
    n=n+1;
% else if we get the julian day not aligned with CTD time point, we will
% need to approximate to closest depth
else
julian=hour/24+minute/24/60+second/24/60/60;
[val,idx]=min(abs(500*(julian_day-floor(julian_day))-500*julian));
idx=idx(1);
if floor(julian_day(end))~=floor(julian_day(1)) && julian>0.8
    julian = julian-1;
end
if julian<julian_day(1)-floor(julian_day(end)) || julian>julian_day(end)-floor(julian_day(end))
   ;
else
    
        Comb{n,1}=table2cell(data_particle(j,1));
        Comb{n,2}=table2cell(data_particle(j,2));
        Comb{n,3}=table2cell(data_particle(j,3));
        Comb{n,4}=table2cell(data_particle(j,4));
        Comb{n,5}=table2cell(data_particle(j,5));
        Comb{n,6}=table2cell(data_particle(j,6));
        Comb{n,7}=table2cell(data_particle(j,7));
        Comb{n,8}=table2cell(data_particle(j,8));
        Comb{n,9}=table2cell(data_particle(j,9));
        Comb{n,10}=table2cell(data_particle(j,10));
        Comb{n,11}=table2cell(data_particle(j,11));
        Comb{n,12}=table2cell(data_particle(j,12));
        Comb{n,13}=table2cell(data_particle(j,13));
        Comb{n,14}=table2cell(data_particle(j,14));
        Comb{n,15}=(data_ctd(idx,1));
        Comb{n,16}=(data_ctd(idx,2));
        Comb{n,17}=(data_ctd(idx,3));
        Comb{n,18}=(data_ctd(idx,4));
        Comb{n,19}=(data_ctd(idx,5));
        Comb{n,20}=(data_ctd(idx,6));
        Comb{n,21}=(data_ctd(idx,7));
        Comb{n,22}=(data_ctd(idx,8));
        Comb{n,23}=(data_ctd(idx,9));
        Comb{n,24}=(data_ctd(idx,10));
        Comb{n,25}=(data_ctd(idx,11));
        Comb{n,26}=(data_ctd(idx,12));
        Comb{n,27}=(data_ctd(idx,13));
        Comb{n,28}=(data_ctd(idx,14));
        Comb{n,29}=(data_ctd(idx,15));
      

        
n=n+1;
end
end

% if it goes to next loop, already write data at the same time point and the next one does not match, then end the loop and for the next time point, start from this row   
end
% pp=vertcat(Comb{:,1};Comb{:,2});
L=cell2table(Comb);
save_path='./RedSea_analysis_test/comb/SPACE/T=29/sub500/cast-5/';
if ~exist(save_path)
  mkdir(save_path)
end
writetable(splitvars(L),[save_path,save_name]);

        
        
toc          

        
        
        
        
        