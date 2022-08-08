%%
    clear all;
    close all;

    myFolder = ('F:\diff\cast-2\');
    if ~isdir(myFolder)
    errorMessage = sprintf('Error: The following folder does not exist:%s', myFolder);
    uiwait(warndlg(errorMessage));
    return;
    end
    path='F:\seg\T13\cast-2\';
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
for i = 50911: length(FileList)
    
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
    imwrite(ROI, [save_path,'/cast-2_image', num2str(image_number),'_',num2str(k), '_', dat,'_',tim,'_', loc,'.tiff'])
%     end
end

% zip(save_path,save_path);
% rmdir(save_path,'s');
end


        