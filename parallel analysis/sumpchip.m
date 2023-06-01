

%%
% interpolate the missing values in some depths

clear
path='/home/hhuan006/sum_2021fall/mi>20/T=1/';
list=dir(fullfile(path,'*cast-1*'));
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


data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

if min(ind_valid)~=1
data_particle(1:min(ind_valid),2)=array2table(nan([0:min(ind_valid)-1,1]));
data_particle(1:min(ind_valid),3)=array2table(nan([0:min(ind_valid)-1,1]));
data_particle(1:min(ind_valid),4)=array2table(nan([0:min(ind_valid)-1,1]));
end

% if max(ind_valid)~=full_row
% data_particle(max(ind_valid)+1:end,2)=array2table(nan([full_row-max(ind_valid),1]));
% data_particle(max(ind_valid)+1:end,3)=array2table(nan([full_row-max(ind_valid),1]));
% data_particle(max(ind_valid)+1:end,4)=array2table(nan([full_row-max(ind_valid),1]));
% end

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
 save_path='/home/hhuan006/sum_2021fall/pchip/mi>20/T=1/';
if ~exist(save_path)
  mkdir(save_path)
end

data_particle(:,2)=array2table(number(:));
data_particle(:,3)=array2table(avevol(:));
data_particle(:,4)=array2table(totvol(:));

% if max(ind_valid)~=full_row
% data_particle(max(ind_valid)+1:end,2)=array2table(nan([full_row-max(ind_valid),1]));
% data_particle(max(ind_valid)+1:end,3)=array2table(nan([full_row-max(ind_valid),1]));
% data_particle(max(ind_valid)+1:end,4)=array2table(nan([full_row-max(ind_valid),1]));
% end

data_particle_up=data_particle;

data_particle=[data_particle_down;data_particle_up];
end


writetable((data_particle),[save_path 'sumpchip_' a{2} '_' a{3}]);
end

        

       

        
       

        