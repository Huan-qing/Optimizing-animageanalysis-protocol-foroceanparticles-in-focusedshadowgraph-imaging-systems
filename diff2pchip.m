clear all;close all


% save the path of images to a variable
file_path =  '/home/hhuan006/2021camera/';
save_dir = '/home/hhuan006/diff2021/';
% make all the images into a list
FileList   = dir(fullfile(file_path, '*cast-7*'));
% sort the images in folder as an ascending order(the order is not right for the initial order)
[~, Index] = natsort({FileList.name});
% return the images in new order into the list
FileList   = FileList(Index);
% gain the entire number of the image files
for ind=1%:length(FileList)
    cast_name=FileList(ind).name;
    subdir=strcat(file_path,cast_name);
    subList   = dir(fullfile(subdir, '*-*'));
% sort the images in folder as an ascending order(the order is not right for the initial order)
    [~, Index] = natsort({subList.name});
    % return the images in new order into the list
    subList   = subList(Index);

L={};
len=length(subList);
image_name=subList(1).name;
I=im2double((imread(strcat(subdir,'/',image_name))));
blankImage=zeros(size(I));
% %imshow(I2_gray)
img_num=100;
for i=len-99:len  % change 5 to 100 later
    image_name = subList(i).name;
    blank=im2double((imread(strcat(subdir,'/',image_name))));
    % First do a check for matching rows, columns, and number of color channels.  Then:
    blankImage = blankImage +  blank/img_num;
end
%------load the artifact-------------------------
    f=1;
save_path = strcat(save_dir,cast_name);
if ~exist(save_path)
    mkdir(save_path)
end

% cd '/home/hhuan006/Bermuda_camera/Bermuda_FoSi2_CTD_camera_April_2019/OFPCTD1/'

for j = 1:2:length(subList)-1

    % create folder to save images
%     if ~exist(path)
%     % creat a new file if not exist
%     mkdir(path)
%     end
    p=0;
    % num_thre=zeros(length(60:3:0),2);
    % get the image name of the 2nd image in the list
    image_name = subList(j).name;
    %split the name based on ' '
    a = regexp(image_name,' ','split');
    % get the variable of image number
    image_label = a{1};
    % read the number after image
    image_number=sscanf(image_label,'Image%f');
    % get the image name of the 2nd image in subtraction
    new=image_name(end-21:end-5);
    % filter the image with even number   image_name_ms = FileList_match(1).name;
    image_name_1 = subList(j).name;
    I1=im2double(imread(strcat(subdir,'/',image_name_1)));
    %imshow(I1)
    %print('-dbmp',['D:/DQ/DQ/MS/research/image analysis/random images/2019spr/random10/canny+sobel+sim/ori_cast-2_1'])
    image_name_2 = subList(j+1).name;
    I2=im2double(imread(strcat(subdir,'/',image_name_2)));
    %imshow(I2)
    %print('-dbmp',['D:/DQ/DQ/MS/research/image analysis/random images/2019spr/random10/canny+sobel+sim/ori_cast-2_2'])

    [cor_I1]=correction_function_noclip(blankImage, I1);
    [cor_I2]=correction_function_noclip(blankImage, I2);
    % %%%figure
    % %imshow(cor_I1)
    % %%%figure
    % %imshow(cor_I2)
    diff=imabsdiff(cor_I1,cor_I2);
    % 
    %%%figure
%     %diff(art_mask==1)=0;
    %imshow(diff)
    imwrite(diff,[save_path '/diff_gain' num2str(image_number) ' ' num2str(new) '.tiff']);
    
% name=['D:/DQ/DQ/MS/research/image analysis/random images/2019spr/new canny/threshold_step/cast-7/image' num2str(image_number) ' ' num2str(new) '.mat']
% save(name,'num_thre')
% clear num_thre
f=f+1;
end
end



%%
    clear all;
    close all;

    myFolder = ('/home/hhuan006/diff2021/cast-7/');
    if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
    end
    path='/home/hhuan006/seg2021/T13/cast-7_T13/';
    if ~exist( path)
    % creat a new file if not exist
    mkdir( path)
    end
    % Change to whatever pattern you need.
    FileList   = dir(fullfile(myFolder,'*.tif*'));
    % sort the images in folder as an ascending order(the order is not right for the initial order)
    [~, Index] = natsort({FileList.name});
    % return the images in new order into the list
    FileList   = FileList(Index);
%     art_mask = load('2019_all_artifact_0902.mat');
%     art_mask = art_mask.close_art;

    n=0;
    M={};
%     cd '/home/hhuan006/Bermuda_subtracted_images/cast-492 subtracted'
for i = 1: length(FileList)
    baseFileName = FileList(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    %fprintf(1, 'Now reading %s', fullFileName);
    
    I = im2double(imread(fullFileName));
%     I(art_mask==1)=0;
    % set the threshold for binary
    brightness_threshold = 13;
    % get the path, name(mainly we get) and file extensions
    [filepath,name,ext] = fileparts(baseFileName);
    % seperate the data based on the ' ' into image number, date and time
    a = regexp(name,' ','split');
    % get the variable of image number
    image_label = a{1};
    % read the number after image
    [image_number,~,~]=sscanf(image_label,'diff_gain%f');
    % read the date on the name of image
    % read the date on the name of image
%     dat=baseFileName(end-16:end-9);
    dat=baseFileName(end-21:end-14);
    % read the time on the name of image
%     tim=baseFileName(end-7:end);
    tim=baseFileName(end-12:end-5);
    % change the '-' to ':' on name of time
    tim= strrep(tim,':','-');
    new=baseFileName(end-16:end);
    folder = fileparts(which(baseFileName)); % Determine where demo folder is (works with all versions).
    fullFileName = fullfile(folder, baseFileName);
  
    % Threshold the image to get a binary image (only 0's and 1's) of class "logical."
    % Method %1: using im2bw()
    %   normalizedThresholdValue = 0.4; % In range 0 to 1.
    %   thresholdValue = normalizedThresholdValue * max(max(originalImage)); % Gray Levels.
    %   binaryImage = im2bw(originalImage, normalizedThresholdValue);       % One way to threshold to binary
    % Method %2: using a logical operation.
    thresholdValue = brightness_threshold;
    binaryImage = rem_simcan(I, thresholdValue/255); % Bright objects will be chosen if you use >.

%     %imshow(binaryImage)
% ========== IMPORTANT OPTION ============================================================
% Use < if you want to find dark objects instead of bright objects.
%   binaryImage = originalImage < thresholdValue; % Dark objects will be chosen if you use <.
% Do a "hole fill" to get rid of any background pixels or "holes" inside the blobs.
% Get all the blob properties.  Can only pass in originalImage in version R2008a and later.
blobMeasurements = regionprops(binaryImage, I, 'all');%ABB
numberOfBlobs = size(blobMeasurements, 1);
% textFontSize = 14;    % Used to control size of "blob number" labels put atop the image.
% labelShiftX = -7; % Used to align the labels in the centers of the coins.
% blobECD = zeros(1, numberOfBlobs);
% Print header line in the command window.
% %fprintf(1,'Image% Blob%      Mean Intensity  Area   Perimeter    Centroid       Diameter');

save_path= [path '/image-' num2str(image_number)];
if ~exist(save_path)
    mkdir(save_path);
end
%M=cell(numberOfBlobs,14); %preload matrix M with zeros
% Loop over all blobs printing their measurements to the command window.
for k = 1 : numberOfBlobs           % Loop through all blobs.
%     area=blobMeasurements(k).Area;
%     if area>50
    location = blobMeasurements(k).BoundingBox;
    loc=['loc_',num2str(location(1)),'_',num2str(location(2)),'_',num2str(location(3)),'_',num2str(location(4))];
    ROI=imcrop(I, location);
%     ROI=adapthisteq(ROI);
%     %imshow(ROI)
%     location=
    imwrite(ROI, [save_path,'/cast-7_image', num2str(image_number),'_',num2str(k), '_', dat,'_',tim,'_', loc,'.tiff'])
%     end
end
zip(save_path,save_path);
rmdir(save_path,'s');
end


% date_time data
clear all;
    close all;
    cd '/home/hhuan006/diff2021/cast-5/'
    myFolder = ('/home/hhuan006/diff2021/cast-5/');
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
for i = 1: length(theFiles)
    baseFileName = theFiles(i).name;
    fullFileName = fullfile(myFolder, baseFileName);
    fprintf(1, 'Now reading %s', fullFileName);
    

    % set the threshold for binary
    brightness_threshold = 13;
    % get the path, name(mainly we get) and file extensions
    [filepath,name,ext] = fileparts(baseFileName);
    % seperate the data based on the ' ' into image number, date and time
    a = regexp(name,' ','split');
    % get the variable of image number
    image_label = a{1};
    % read the number after image
    [image_number,~,~]=sscanf(image_label,'diff_gain%f');
    % read the date on the name of image
    % read the date on the name of image
    dat_{i}=baseFileName(end-21:end-13);
    % read the time on the name of image
    tim_{i}=baseFileName(end-12:end-5)
    % change the '-' to ':' on name of time
    tim_{i}= strrep(tim_{i},'-',':');
    new=baseFileName(end-21:end);
    folder = fileparts(which(baseFileName)); % Determine where demo folder is (works with all versions).
    fullFileName = fullfile(folder, baseFileName);
    originalImage = imread(fullFileName);
    % Check to make sure that it is grayscale, just in case the user substituted their own image.
    [rows, columns, numberOfColorChannels] = size(originalImage);
    thresholdValue = brightness_threshold;
    binaryImage = rem_simcan(originalImage,  thresholdValue); 
    imshow(binaryImage)

blobMeasurements = regionprops(binaryImage, originalImage, 'all');%ABB
numberOfBlobs = size(blobMeasurements, 1);

fprintf(1,'Image% Blob%      Mean Intensity  Area   Perimeter    Centroid       Diameter');
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
savepath='/home/hhuan006/batch2021spr/T13/';
if ~exist(savepath)
  mkdir(savepath)
end
writetable(L,[savepath 'T13_date_time-cast-5.csv']); 


%%


% combine with CTD data

clear; close all; clc
% go to the path where the files is, and load the csv
path= '/home/hhuan006/batch2021spr/T13/';
path_ctd='/home/hhuan006/2021_Spr_CTD/';
ctd_list=dir(fullfile(path_ctd,'*cast-5*'));
ctd_name=ctd_list(1).name;
data_ctd=csvread(strcat(path_ctd,ctd_name),1,0);
filelist=dir(fullfile(path,'*cast-5.csv'));
n=1;
Comb={};
m=1;
julian_day=data_ctd(:,17);
% get the row number and column number of the CTD data
[hour_ctd,minute_ctd,second_ctd]=julian_conversion(julian_day);
L={};
% L = num2cell(zeros(0,34));
% L.Properties.VariableNames = {'image_number', 'k', 'Area', 'MajorAxisLength',  'MinorAxisLength', 'EquivDiameter', 'Eccentricity', 'Orientation', 'FilledArea', 'Solidity', 'Perimeter', 'MeanIntensity','date','time', 'Conductivity','Conductivity2','temperature', 'temperature2','pressure', 'oxy', 'oxy2','fluorescence[ug/l]',  'fluorescence[mg/m^3]', 'BeamAttenuation','BeamTransmission', 'PSU', 'PSU2','density','density2','depth','days','hours','min','sec'}
file_name=filelist(1).name;    
data_particle=readtable(strcat(path,file_name));
% get the time in image data
time_point=table2cell(data_particle(:,14));
[w,p]=size(time_point);
% change the format of x into string as B
time_p=string(time_point);
for j=1:length(time_point)
    tic



% give the data of Julian days into a variable A

% calculate the julian days into hour, minute and second, and save into
% seperate variables

% split the string in the image data into hour, minute and second
unt=time_p(j,1);
P=strsplit(unt,':');
hour=str2num(P{1});
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
        Comb{n,15}=num2cell(data_ctd(i,1));
        Comb{n,16}=num2cell(data_ctd(i,2));
        Comb{n,17}=num2cell(data_ctd(i,3));
        Comb{n,18}=num2cell(data_ctd(i,4));
        Comb{n,19}=num2cell(data_ctd(i,5));
        Comb{n,20}=num2cell(data_ctd(i,6));
        Comb{n,21}=num2cell(data_ctd(i,7));
        Comb{n,22}=num2cell(data_ctd(i,8));
        Comb{n,23}=num2cell(data_ctd(i,9));
        Comb{n,24}=num2cell(data_ctd(i,10));
        Comb{n,25}=num2cell(data_ctd(i,11));
        Comb{n,26}=num2cell(data_ctd(i,12));
        Comb{n,27}=num2cell(data_ctd(i,13));
        Comb{n,28}=num2cell(data_ctd(i,14));
        Comb{n,29}=num2cell(data_ctd(i,15));
        Comb{n,30}=num2cell(data_ctd(i,16));
        Comb{n,31}=num2cell(data_ctd(i,17));
        Comb{n,32}=num2cell(data_ctd(i,18));
        Comb{n,33}=num2cell(data_ctd(i,19));
        Comb{n,34}=num2cell(data_ctd(i,20));
    
  
        %[image_number k Area MajorAxisLength  MinorAxisLength EquivDiameter Eccentricity Orientation FilledArea Solidity Perimeter MeanIntensity data time Conductivity   Conductivity2  temperature temperature2   pressure    oxy oxy2   fluorescence[ug/l]  fluorescence[mg/m^3]  BeamAttenuation  BeamTransmission PSU PSU2   density density2   depth   days    hours   min  sec ]
    n=n+1;
% else if we get the julian day not aligned with CTD time point, we will
% need to approximate to closest depth
else
julian=hour/24+minute/24/60+second/24/60/60;
[val,idx]=min(abs(100*(julian_day-floor(julian_day))-100*julian));
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
        Comb{n,15}=num2cell(data_ctd(idx,1));
        Comb{n,16}=num2cell(data_ctd(idx,2));
        Comb{n,17}=num2cell(data_ctd(idx,3));
        Comb{n,18}=num2cell(data_ctd(idx,4));
        Comb{n,19}=num2cell(data_ctd(idx,5));
        Comb{n,20}=num2cell(data_ctd(idx,6));
        Comb{n,21}=num2cell(data_ctd(idx,7));
        Comb{n,22}=num2cell(data_ctd(idx,8));
        Comb{n,23}=num2cell(data_ctd(idx,9));
        Comb{n,24}=num2cell(data_ctd(idx,10));
        Comb{n,25}=num2cell(data_ctd(idx,11));
        Comb{n,26}=num2cell(data_ctd(idx,12));
        Comb{n,27}=num2cell(data_ctd(idx,13));
        Comb{n,28}=num2cell(data_ctd(idx,14));
        Comb{n,29}=num2cell(data_ctd(idx,15));
        Comb{n,30}=num2cell(data_ctd(idx,16));
        Comb{n,31}=num2cell(data_ctd(idx,17));
        Comb{n,32}=num2cell(data_ctd(idx,18));
        Comb{n,33}=num2cell(data_ctd(idx,19));
        Comb{n,34}=num2cell(data_ctd(idx,20));
        
n=n+1;
end
end

% if it goes to next loop, already write data at the same time point and the next one does not match, then end the loop and for the next time point, start from this row   
toc
end
% pp=vertcat(Comb{:,1};Comb{:,2});
L=cell2table(Comb);
save_path='/home/hhuan006/comb2021spr/T13/';
if ~exist(save_path)
  mkdir(save_path)
end
writetable(splitvars(L),[save_path,'combine_cast-5_T13.csv']);

%%
% merge into depth wise data in 1m

  %image_number k Area MajorAxisLength  MinorAxisLength EquivDiameter Eccentricity Orientation FilledArea Solidity Perimeter MeanIntensity Conductivity   Conductivity2  temperature temperature2   pressure    oxy oxy2   fluorescence[ug/l]  fluorescence[mg/m^3]  BeamAttenuation  BeamTransmission PSU PSU2   density density2   depth   days    hours   min  sec 
clear all; close all; clc
path='/home/hhuan006/comb2021spr/T13/';
list=dir(fullfile(path,'*cast-5*'));
for p=1:length(list)
    file_name=list(p).name;
    a = regexp(file_name,'_','split');
    thre_num=a{3};
    thre_num=thre_num(1:end-4);
    data=readtable(strcat(path,file_name));
    dat=table2cell(data);
    [cc,rr]=size(dat);
% filter the image with particle area>1000, record the row number of Total
% Volume>1000, and delete all.
area=table2array(data(:,3));
dat=dat(area<=1000,:);
%get the size of matrix again and put features into different variables,
%change area back to radius and calculate the volume again
[cc,rr]=size(dat);
num=dat(:,1);
k=dat(:,2);
area=cell2mat(dat(:,3));
volume=4/3.*pi.*((area/pi).^0.5).^3;
 
par_ctd=dat(:,4:34);
 
n=1;
image_aver={};
number=1;
mj=0;
num_par=1;
num_im=1;
totvol=0;
pp=1;
for i= 2:cc %cc 
 %calculate and average the particle number, average volume and total volume of the same
 %depth to one image to represent the data of the depth
 if cell2mat(par_ctd(i-1,27))~=cell2mat(par_ctd(i,27))
     
      if cell2mat(par_ctd(i-1,27))==cell2mat(par_ctd(i-2,27))
        if cell2mat(num(i-2,1))==cell2mat(num(i-1,1))
          
           totvol=totvol+volume(i-1);  %to get the sum volume per image
           mj=totvol/num_par+mj;
           num_par=1;
           totvol=0;
        else 
           totvol=totvol+volume(i-1);
           mj=totvol+mj;
           num_par=1;
           totvol=0;
           num_im=num_im+1;
        end
      end
  
     
           
       [image_aver{n,1}]=par_ctd(i-1,27);
       [image_aver{n,2}]=num2cell(number/num_im); %number per image
       [image_aver{n,3}]=num2cell(mj/num_im);% calculate the volume
       [image_aver{n,4}]=num2cell(number/num_im*mj/num_im);% calculate the total volume
       [image_aver{n,5}]=par_ctd(i-1,10);
       [image_aver{n,6}]=par_ctd(i-1,11);
       [image_aver{n,7}]=par_ctd(i-1,12);
       [image_aver{n,8}]=par_ctd(i-1,13);
       [image_aver{n,9}]=par_ctd(i-1,14);
       [image_aver{n,10}]=par_ctd(i-1,15);
       [image_aver{n,11}]=par_ctd(i-1,16);
       [image_aver{n,12}]=par_ctd(i-1,17);
       [image_aver{n,13}]=par_ctd(i-1,18);
       [image_aver{n,14}]=par_ctd(i-1,19);
       [image_aver{n,15}]=par_ctd(i-1,20);
       [image_aver{n,16}]=par_ctd(i-1,21);
       [image_aver{n,17}]=par_ctd(i-1,22);
       [image_aver{n,18}]=par_ctd(i-1,23);
       [image_aver{n,19}]=par_ctd(i-1,24);
       [image_aver{n,20}]=par_ctd(i-1,25);
       [image_aver{n,21}]=par_ctd(i-1,26);
       [image_aver{n,22}]=par_ctd(i-1,28);
       [image_aver{n,23}]=par_ctd(i-1,29);
       [image_aver{n,24}]=par_ctd(i-1,30);
       [image_aver{n,25}]=par_ctd(i-1,31);
    

 
       number=1;
       n=n+1;
       mj=0;
       num_im=1; 
  % if there is an image with 0 particle then go to the next cycle
  elseif cell2mat(k(i))==0  
          continue
   else
       number=number+1;
%        mj=mj+volume(i-1);
       if cell2mat(num(i,1))==cell2mat(num(i-1,1))
           num_par=num_par+1;
           totvol=totvol+volume(i-1);  %get the average volume per image
       else 
           totvol=totvol+volume(i-1);
           mj=totvol/num_par+mj;
           num_par=1;
           totvol=0;
           num_im=num_im+1;
       end
   end
 
end
 
 
% to calculate the last time point, because in cycle i, we input the data
% of cycle i-1
       [image_aver{n,1}]=par_ctd(i-1,27);
       [image_aver{n,2}]=num2cell(number/num_im); %number per image
       [image_aver{n,3}]=num2cell(mj/num_im);% calculate the volume
       [image_aver{n,4}]=num2cell(number/num_im*mj/num_im);% calculate the total volume
       [image_aver{n,5}]=par_ctd(i-1,10);
       [image_aver{n,6}]=par_ctd(i-1,11);
       [image_aver{n,7}]=par_ctd(i-1,12);
       [image_aver{n,8}]=par_ctd(i-1,13);
       [image_aver{n,9}]=par_ctd(i-1,14);
       [image_aver{n,10}]=par_ctd(i-1,15);
       [image_aver{n,11}]=par_ctd(i-1,16);
       [image_aver{n,12}]=par_ctd(i-1,17);
       [image_aver{n,13}]=par_ctd(i-1,18);
       [image_aver{n,14}]=par_ctd(i-1,19);
       [image_aver{n,15}]=par_ctd(i-1,20);
       [image_aver{n,16}]=par_ctd(i-1,21);
       [image_aver{n,17}]=par_ctd(i-1,22);
       [image_aver{n,18}]=par_ctd(i-1,23);
       [image_aver{n,19}]=par_ctd(i-1,24);
       [image_aver{n,20}]=par_ctd(i-1,25);
       [image_aver{n,21}]=par_ctd(i-1,26);
       [image_aver{n,22}]=par_ctd(i-1,28);
       [image_aver{n,23}]=par_ctd(i-1,29);
       [image_aver{n,24}]=par_ctd(i-1,30);
       [image_aver{n,25}]=par_ctd(i-1,31);
% change the format of variable to make the cell all cell
image_aver=cell2table(image_aver);
image_aver=table2cell(image_aver);
% read the CTD data and convert it into cell format
ori_ctd=readtable(['/home/hhuan006/2021_Spr_CTD/' a{2} '.csv']);
julian_ctd=table2array(ori_ctd(:,17));
ori_ctd=table2cell(ori_ctd);
[h,l]=size(ori_ctd);
[hh,ll]=size(image_aver);
sum={};

% list all the data in CTD with depth on the first column
for i=1:h
    
       [sum{i,1}]=ori_ctd(i,16);
       [sum{i,2}]=NaN;
       [sum{i,3}]=NaN;
       [sum{i,4}]=NaN;
       [sum{i,5}]=NaN;
       [sum{i,6}]=NaN;
       [sum{i,7}]=ori_ctd(i,1);
       [sum{i,8}]=ori_ctd(i,2);
       [sum{i,9}]=ori_ctd(i,3);
       [sum{i,10}]=ori_ctd(i,4);
       [sum{i,11}]=ori_ctd(i,5);
       [sum{i,12}]=ori_ctd(i,6);
       [sum{i,13}]=ori_ctd(i,7);
       [sum{i,14}]=ori_ctd(i,8);
       [sum{i,15}]=ori_ctd(i,9);
       [sum{i,16}]=ori_ctd(i,10);
       [sum{i,17}]=ori_ctd(i,11);
       [sum{i,18}]=ori_ctd(i,12);
       [sum{i,19}]=ori_ctd(i,13);
       [sum{i,20}]=ori_ctd(i,14);
       [sum{i,21}]=ori_ctd(i,15);
       [sum{i,22}]=ori_ctd(i,17);
       [sum{i,23}]=ori_ctd(i,18);
       [sum{i,24}]=ori_ctd(i,19);
       [sum{i,25}]=ori_ctd(i,20);
     
end
% add those match in image data to the CTD variable
for i=1:h
    for j=1:hh
      if cell2mat(ori_ctd(i,17))==cell2mat(image_aver(j,22))
       [hour,minute,second]=julian_conversion(julian_ctd(i,:));
       [sum{i,2}]=image_aver(j,2);
       [sum{i,3}]=image_aver(j,3);
       [sum{i,4}]=image_aver(j,4);
       [sum{i,5}]=image_aver(j,5);
       [sum{i,6}]=[num2str(hour),':',num2str(minute),':',num2str(second)];
      
      end
   end
end
save_path='/home/hhuan006/sum2021spr/raw/T13';
if ~exist(save_path)
  mkdir(save_path)
end
L=cell2table(sum);
writetable(L,[save_path '/sum_' a{2} '_' a{3}]);
end

%%
% interpolate the missing values in some depths

clear
path='/home/hhuan006/sum2021spr/raw/T13/';
list=dir(fullfile(path,'*cast-5*'));
m=0;
for p=1:length(list)
    file_name=list(p).name
    a = regexp(file_name,'_','split');
    thre_num=a{3};
    thre_num=strsplit(thre_num,'.');
    thre_num=thre_num{1};
    data=readtable(strcat(path,'/',file_name));
    depth=double(table2array(data(:,1)));
ind_max_depth=find(depth==max(depth));
% pchip for the downward cast first
data_particle=data(1:ind_max_depth,:);
depth=table2array(data_particle(:,1));  % read the data of depth
number=table2array(data_particle(:,2));  % read the data of number
avevol=table2array(data_particle(:,3));  % read the data of avevol
totvol=table2array(data_particle(:,4));  % read the data of totvol
if iscell(depth)
    depth=str2double(depth);
end
if iscell(number)
    number=str2double(number);
end
if iscell(avevol)
    avevol=str2double(avevol);
end
if iscell(totvol)
    totvol=str2double(totvol);
end
ind=find(depth==max(depth));

% depth=cell2mat(data_particle(:,1));  % read the data of depth
% number=cell2mat(data_particle(:,2));  % read the data of number
% avevol=cell2mat(data_particle(:,3));  % read the data of avevol
% totvol=cell2mat(data_particle(:,4));  % read the data of 
% f=cell2mat(data_particle(:,5));

[ss,pp]=size(depth);
 
% if there are any NaN value in any data, use pchip to fill a number based
% on other data
for i=1:ss
    if ~isnan(number(i))
        m=1;
    elseif m==1
  if isnan(number(i))
     number(i)=pchip(depth,number,depth(i));
     number(i);
 
  end
  if isnan(avevol(i))
    avevol(i)=pchip(depth,avevol,depth(i));
    avevol(i);
 
  end
  if isnan(totvol(i))
    totvol(i)=pchip(depth,totvol,depth(i));
    totvol(i);
 
 
  end
    end
    % if there are any data<0, which is not realistic, then leave it as NaN
  if number(i)<0| avevol(i)<0| totvol(i)<0
    number(i)=nan;
    avevol(i)=nan;
    totvol(i)=nan;
    
    end
end
 save_path='/home/hhuan006/sum2021spr/pchip/T13/';
if ~exist(save_path)
  mkdir(save_path)
end
data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

data_particle_down=data_particle;

if ind_max_depth<size(data,1);
% pchip for the upward cast later
data_particle=data(ind_max_depth+1:end,:);
depth=table2array(data_particle(:,1));  % read the data of depth
number=table2array(data_particle(:,2));  % read the data of number
avevol=table2array(data_particle(:,3));  % read the data of avevol
totvol=table2array(data_particle(:,4));  % read the data of totvol
if iscell(depth)
    depth=str2double(depth);
end
if iscell(number)
    number=str2double(number);
end
if iscell(avevol)
    avevol=str2double(avevol);
end
if iscell(totvol)
    totvol=str2double(totvol);
end


[ss,pp]=size(depth);
 
% if there are any NaN value in any data, use pchip to fill a number based
% on other data
for i=1:ss
    if ~isnan(number(i))
        m=1;
    elseif m==1
  if isnan(number(i))
     number(i)=pchip(depth,number,depth(i));
     number(i);
 
  end
  if isnan(avevol(i))
    avevol(i)=pchip(depth,avevol,depth(i));
    avevol(i);
 
  end
  if isnan(totvol(i))
    totvol(i)=pchip(depth,totvol,depth(i));
    totvol(i);
 
 
  end
    end
    % if there are any data<0, which is not realistic, then leave it as NaN
  if number(i)<0| avevol(i)<0| totvol(i)<0
    number(i)=nan;
    avevol(i)=nan;
    totvol(i)=nan;
    
    end
end
 save_path='/home/hhuan006/sum2021spr/pchip/T13/';
if ~exist(save_path)
  mkdir(save_path)
end
data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

data_particle_up=data_particle;

data_particle=[data_particle_down;data_particle_up];
end

writetable((data_particle),[save_path 'sumpchip_' a{2} '_' a{3}]);
end

        

       

        
       