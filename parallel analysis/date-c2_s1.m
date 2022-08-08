

% date_time data
clear all;
    close all;
%     cd '/home/hhuan006/diff2021fall/cast-2/'
    myFolder = ('/home/hhuan006/diff2021fall/cast-2/');
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
    sec_total=100;
    sec=1;
    each=fix(length(theFiles)/sec_total);
     % to 11
    if sec==sec_total
        last_file=length(theFiles);
    else
        last_file=each*sec;
    end
for i=each*(sec-1)+1:last_file
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
    tim_{i}=baseFileName(end-12:end-5);
    % change the '-' to ':' on name of time
    tim_{i}= strrep(tim_{i},'-',':');
    new=baseFileName(end-21:end);
    folder = fileparts(which(baseFileName)); % Determine where demo folder is (works with all versions).
    fullFileName = fullfile(folder, baseFileName);
    originalImage = imread(strcat(myFolder,fullFileName));
    % Check to make sure that it is grayscale, just in case the user substituted their own image.
    [rows, columns, numberOfColorChannels] = size(originalImage);
    thresholdValue = brightness_threshold;
    binaryImage = rem_simcan(originalImage,  thresholdValue); 
%     imshow(binaryImage)

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
savepath='/home/hhuan006/batch2021fall/T13/sub100/';
if ~exist(savepath)
  mkdir(savepath)
end
writetable(L,[savepath 'T13_date_time-cast-2_sec=1.csv']); 

        

        