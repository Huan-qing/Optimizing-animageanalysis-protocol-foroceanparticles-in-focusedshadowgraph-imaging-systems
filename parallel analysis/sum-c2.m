%%
% merge into depth wise data in 1m

  %image_number k Area MajorAxisLength  MinorAxisLength EquivDiameter Eccentricity Orientation FilledArea Solidity Perimeter MeanIntensity Conductivity   Conductivity2  temperature temperature2   pressure    oxy oxy2   fluorescence[ug/l]  fluorescence[mg/m^3]  BeamAttenuation  BeamTransmission PSU PSU2   density density2   depth   days    hours   min  sec 
clear all; close all; clc
path='/home/hhuan006/comb2021fall/T18/merge/';
list=dir(fullfile(path,'*cast-2*'));
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
esd=cell2mat(dat(:,6));
volume=4/3.*pi.*(esd/2).^3;
 
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
ori_ctd=readtable(['/home/hhuan006/CTD2021fall/' a{2} '.csv']);
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
save_path='/home/hhuan006/sum2021fall/raw/T18';
if ~exist(save_path)
  mkdir(save_path)
end
L=cell2table(sum);
writetable(L,[save_path '/sum_' a{2} '_' a{3}]);
end

%%
% interpolate the missing values in some depths

clear
path='/home/hhuan006/sum2021fall/raw/T18/';
list=dir(fullfile(path,'*cast-2*'));
m=0;
for p=1:length(list)
    file_name=list(p).name
    a = regexp(file_name,'_','split');
    thre_num=a{3};
    thre_num=strsplit(thre_num,'.');
    thre_num=thre_num{1};
    data=readtable(strcat(path,'/',file_name));
    [full_row,~]=size(data);
    depth=double(table2array(data(:,1)));
ind_max_depth=find(depth==max(depth));
% pchip for the downward cast first
data_particle=data(1:ind_max_depth,:);
depth=table2array(data_particle(:,1));  % read the data of depth
number=table2array(data_particle(:,2));  % read the data of number
avevol=table2array(data_particle(:,3));  % read the data of avevol
totvol=table2array(data_particle(:,4));  % read the data of totvol
ind_valid=find(~isnan(number));
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
 save_path='/home/hhuan006/sum2021fall/pchip/T18/';
if ~exist(save_path)
  mkdir(save_path)
end

data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

if min(ind_valid)~=1
data_particle(1:min(ind_valid),2)=array2table(nan([0:min(ind_valid)-1,1]));
data_particle(1:min(ind_valid),3)=array2table(nan([0:min(ind_valid)-1,1]));
data_particle(1:min(ind_valid),4)=array2table(nan([0:min(ind_valid)-1,1]));
end

if max(ind_valid)~=full_row
data_particle(max(ind_valid)+1:end,2)=array2table(nan([full_row-max(ind_valid),1]));
data_particle(max(ind_valid)+1:end,3)=array2table(nan([full_row-max(ind_valid),1]));
data_particle(max(ind_valid)+1:end,4)=array2table(nan([full_row-max(ind_valid),1]));
end

data_particle_down=data_particle;

if ind_max_depth<size(data,1);
% pchip for the upward cast later
data_particle=data(ind_max_depth+1:end,:);
depth=table2array(data_particle(:,1));  % read the data of depth
number=table2array(data_particle(:,2));  % read the data of number
avevol=table2array(data_particle(:,3));  % read the data of avevol
totvol=table2array(data_particle(:,4));  % read the data of totvol
ind_valid=find(~isnan(number));
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
 save_path='/home/hhuan006/sum2021fall/pchip/T18/';
if ~exist(save_path)
  mkdir(save_path)
end

data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

if max(ind_valid)~=full_row
data_particle(max(ind_valid)+1:end,2)=array2table(nan([full_row-max(ind_valid),1]));
data_particle(max(ind_valid)+1:end,3)=array2table(nan([full_row-max(ind_valid),1]));
data_particle(max(ind_valid)+1:end,4)=array2table(nan([full_row-max(ind_valid),1]));
end

data_particle_up=data_particle;

data_particle=[data_particle_down;data_particle_up];
end


writetable((data_particle),[save_path 'sumpchip_' a{2} '_' a{3}]);
end

        

       