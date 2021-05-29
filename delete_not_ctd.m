img_path='F:\OFPApr2021 Images all casts\camera images\';
ctd_path='F:\OFPApr2021 Images all casts\CTD\2104_AE CTD data\CTD_Data\2021\csv file\';
im_list=dir(fullfile(img_path,'*-*'));
for i=1:length(im_list)
data_ctd=readtable(strcat(ctd_path,[im_list(i).name,'.csv']));
theFiles=dir(fullfile(strcat(img_path,im_list(i).name),'*-*'));
[~, Index] = natsort({theFiles.name});
% return the images in new order into the list
theFiles   = theFiles(Index);
% give the data of Julian days into a variable A
julian_day=tb2number(data_ctd(:,17));
start_time=julian_day(1);
start_time=start_time-floor(start_time);
end_time=julian_day(end);
end_time=end_time-floor(end_time);
for ii= 1:length(theFiles)
    
    baseFileName = theFiles(ii).name;
    %calculate the time point
    tim_{ii}=baseFileName(end-12:end-5);
    % change the '-' to ':' on name of time
    tim_{ii}= strrep(tim_{ii},'-',':');

    unt=tim_{ii};
    P=strsplit(unt,':');
    hour(ii,1)=str2num(P{1});
    minute(ii,1)=str2num(P{2});
    second(ii,1)=str2num(P{3});
end    
julian=hour/24+minute/24/60+second/24/60/60;

[~,i_before]=min(abs(100*start_time-100*(julian-floor(julian))))
[~,i_after]=min(abs(100*end_time-100*(julian-floor(julian))))

% if the time points in 2 variabls match, the write down the data
% if hour==hour_ctd_0_down(1)&minute==minute_ctd_0_down(1)&second==second_ctd_0_down(1)
%   i_0_down=1;
% end
a=[1:i_before,i_after:length(theFiles)];
    for j = 1:length(a)
        delete(strcat(img_path,im_list(i).name,'\',theFiles(a(j)).name))
%         copyfile(strcat(img_path,im_list(i).name,'\',theFiles(a(j)).name),save_path_sub)
    end
end