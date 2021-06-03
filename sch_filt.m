%% for whole cast
file_path =  'F:\OFPApr2021 Images all casts\camera images\cast-1\';
save_path = 'F:\OFPApr2021 Images all casts\sch img\cast-1\';
if ~exist(save_path)
    mkdir(save_path)
end
FileList   = dir(fullfile(file_path, '*.tiff'));
% sort the images in folder as an ascending order(the order is not right for the initial order)
[~, Index] = natsort({FileList.name});
% return the images in new order into the list
FileList   = FileList(Index);
% gain the entire number of the image files
L={};
len=length(FileList);
image_name=FileList(1).name;
I=im2double((imread(strcat(file_path,image_name))));
blankImage=zeros(size(I));
% imshow(I2_gray)
img_num=100;
n=0;
for i=len-99:len  % change 5 to 100 later
image_name = FileList(i).name;
blank=im2double((imread(strcat(file_path,image_name))));
% First do a check for matching rows, columns, and number of color channels.  Then:
blankImage = blankImage +  blank/img_num;
end
for j = 1: length(FileList)
%     tic
image_name = FileList(j).name;
I_sch=im2double((imread(strcat(file_path,image_name))));
% subtract the target image - blankimage, keep minus or plus sign
sub = imsubtract(I_sch,blankImage);
% if the target pixel is brighter than our (schlieren is due to light distortion, some regions has more light from refraction)
mask=sub>0.01;
mask=rem_simcan(sub,0.01);
% remove objects with area<30
large=bwareaopen(mask, 30);
label=bwlabel(large);
major = regionprops(large, 'MajorAxisLength');
for i =1 :length(major)
long(i,1)=major(i).MajorAxisLength;
end
minor = regionprops(large, 'MinorAxisLength');
for i =1 :length(minor)
short(i,1)=minor(i).MinorAxisLength;
end
ratio=short./long;
a = regionprops(large, 'Area');
for i =1 :length(a)
area(i,1)=a(i).Area;
end
% create a mask for schlieren
sch=zeros(size(I_sch));
for i = 1: max(label(:))
    if (ratio(i)<0.2&&area(i)>100)||  area(i)>40000
        sch(label==i)=1;
    end
end
label_sch=bwlabel(sch);
if max(label_sch(:))>=5 || max(area)>40000
    copyfile(strcat(file_path,image_name),save_path);
    n=n+1
end
% figure
% imshow(I_sch)
% % figure
% % imshow(mask)
% figure
% imshow(large)
% figure
% imshow(logical(sch))
% toc
end

%% try single
file_path =  'F:\OFPApr2021 Images all casts\camera images\cast-1\';
I_sch=im2double(imread('D:\DQ\DQ\MS\research\optical experiment\schlieren\exp\Image104773 21-04-14 23-55-28.tiff'));
% make all the images into a list

FileList   = dir(fullfile(file_path, '*.tiff'));
% sort the images in folder as an ascending order(the order is not right for the initial order)
[~, Index] = natsort({FileList.name});
% return the images in new order into the list
FileList   = FileList(Index);
% gain the entire number of the image files
L={};
len=length(FileList);
image_name=FileList(1).name;
I=im2double((imread(strcat(file_path,image_name))));
blankImage=zeros(size(I));
% imshow(I2_gray)
img_num=100;
for i=len-99:len  % change 5 to 100 later
image_name = FileList(i).name;
blank=im2double((imread(strcat(file_path,image_name))));
% First do a check for matching rows, columns, and number of color channels.  Then:
blankImage = blankImage +  blank/img_num;
end
% subtract the target image - blankimage, keep minus or plus sign
sub = imsubtract(I_sch,blankImage);
% if the target pixel is brighter than our (schlieren is due to light distortion, some regions has more light from refraction)
mask=sub>0.01;
% mask=rem_simcan(sub,0.01);
% remove objects with area<30
large=bwareaopen(mask, 30);
label=bwlabel(large);
major = regionprops(large, 'MajorAxisLength');
for i =1 :length(major)
long(i,1)=major(i).MajorAxisLength;
end
minor = regionprops(large, 'MinorAxisLength');
for i =1 :length(minor)
short(i,1)=minor(i).MinorAxisLength;
end
ratio=short./long;
a = regionprops(large, 'Area');
for i =1 :length(a)
area(i,1)=a(i).Area;
end
% create a mask for schlieren
sch=zeros(size(I_sch));
for i = 1: max(label(:))
    if (ratio(i)<0.2&&area(i)>300) ||  area(i)>40000
        sch(label==i)=1;
    end
end
label_sch=bwlabel(sch);
if max(label_sch(:))>=5 || max(area)>40000
    'is schlieren'
end
figure
imshow(I_sch)
% figure
% imshow(mask)
figure
imshow(large)
figure
imshow(logical(sch))
% tem=bwareaopen(large, 300);
% figure;imshow(tem)