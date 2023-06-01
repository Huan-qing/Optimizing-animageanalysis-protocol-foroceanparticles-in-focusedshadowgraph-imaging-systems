



tic
clear; close all; clc
path='/home/hhuan006/comb_2021fall/SPACE/T=29/sub100/cast-1/';
% path='/home/hhuan006/comb_2021fall/SPACE/T=29/sub100/cast-1/';
list=dir(fullfile(path,'*cast*'));
[~,index]=natsort({list.name});
list=list(index);
cat=cell(1,length(list));
m_big={};
for t=1:length(list)
    
        file_name=list(t).name;
        a = regexp(file_name,'_','split');
        cs = regexp(a{2},'ne-','split');
        cat=readtable(strcat(path,file_name));
        if height(cat)==0
            continue
        end
        esd=cat{:,6};
        % filter out the particles greater than 1 mm/1000 um
        id_small=find(esd*11.6166<=1000);
        cat=cat(id_small,:);
        img_num=cat{:,1};
        image_idx=findgroups(img_num);
        img_num_uni=unique(img_num);
        particle_idx=cat{:,2};
        particle_no=groupsummary(particle_idx,image_idx,'max');
        volume=4/3.*pi.*(cat{:,6}/2).^3;
        volume_aver=groupsummary(volume,image_idx,'mean');
        total_volume=groupsummary(volume,image_idx,'sum');
        ctd_rest=cat{:,15:34};
        ctd_rest=groupsummary(ctd_rest,image_idx,'mean');

        idx=max(image_idx);
        date={};
        time={};
    for i=1:idx
        id_uni=find(image_idx==i);
        date=[date;cat(id_uni(1),13)];
        time=[time;cat(id_uni(1),14)]; 
    end
        m=[table(img_num_uni,particle_no,volume_aver,total_volume,table2cell(date),table2cell(time)),array2table(ctd_rest)];
        m_big=[m_big;m];
%     if size(cat,1)~=0
%         m=[m;cat];
%     end
end

%% down-cast (duplicate depth in two-way cast)

depth_bin=m_big{:,22};
idx_down=max(find(depth_bin==max(depth_bin)));
m_big_down=m_big(1:idx_down,:);
depth_bin=m_big_down{:,22};
depth_idx=findgroups(depth_bin);
depth_group_uni=unique(depth_idx);
depth_bin_uni=unique(depth_bin);
particle_no_bin=groupsummary(m_big_down{:,2},depth_idx,'mean');
volume_aver_bin=groupsummary(m_big_down{:,3},depth_idx,'mean');
total_volume_bin=groupsummary(m_big_down{:,4},depth_idx,'mean');
ctd_rest_perim=m_big_down{:,7:26};
ctd_rest_bin=groupsummary(ctd_rest_perim,depth_idx,'mean');
ctd_rest_bin(:,16)=[];
date_bin={};
time_bin={};

for i = depth_group_uni'
        id_bin=find(depth_idx==i);
        date_bin=[date_bin;m_big_down(id_bin(1),5)];
        time_bin=[time_bin;m_big_down(id_bin(1),6)]; 
end
    
m_bin_ori_down=[table(depth_bin_uni,particle_no_bin,volume_aver_bin,total_volume_bin,table2cell(date_bin),table2cell(time_bin)),array2table(ctd_rest_bin)];
name_col={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time','Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2','days','hours','min','sec'};
m_bin_ori_down=renamevars(m_bin_ori_down,1:width(m_bin_ori_down),name_col);

%%
% find missing data in the table
ori_ctd=csvread(['/home/hhuan006/CTD2021fall/' cs{2} '.csv'],1,0);
% ori_ctd=readtable(['/home/hhuan006/CTD2021fall/' cs{2} '.csv'],1,0);

miss_dep=[];
m_bin_down=m_bin_ori_down;
dp=m_bin_down{:,1};
dp_ctd=ori_ctd(:,16);
idx_max_dp_ctd=find(dp_ctd==max(dp));
dp_ctd_down=dp_ctd(1:idx_max_dp_ctd);
name_col={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time','Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2','days','hours','min','sec'};
name_col_1={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time'};
name_col_2={'Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2'};
name_col_3={'days','hours','min','sec'};

if dp(length(dp))-dp(1)+1~=length(dp)
    for i = 1:length(dp_ctd_down)-1
        dp=m_bin_down{:,1};
        if dp(i)~=dp(i+1)-1
            
            miss_dep=[miss_dep;dp(i)+1];
            
            point_ms=dp(i)+1;
            id_ms_ctd=find(dp_ctd_down==point_ms);
            id_ms_img=min(find(dp==dp(i)));
            pre=m_bin_down(1:id_ms_img,:);
            post=m_bin_down(id_ms_img+1:height(m_bin_down),:);
            t=renamevars(array2table(ori_ctd(id_ms_ctd,17:20)),1:4,name_col_3);
            column_ms=[array2table([point_ms,NaN,NaN,NaN,NaN,NaN],'VariableNames',name_col_1),array2table(ori_ctd(id_ms_ctd,1:15),'VariableNames',name_col_2),t];

            m_bin_down=[pre;column_ms;post];
        end
    end
end

if height(m_big)==height(m_big_down)
    m_bin=[m_bin_down];
    return;
else
%% up-cast (duplicate depth in two-way cast)

depth_bin=m_big{:,22};
idx_up=max(find(depth_bin==max(depth_bin)));
m_big_up=m_big(idx_up+1:length(depth_bin),:);
depth_bin=m_big_up{:,22};
depth_idx=findgroups(depth_bin);
depth_group_uni=unique(depth_idx);
depth_bin_uni=unique(depth_bin);
particle_no_bin=groupsummary(m_big_up{:,2},depth_idx,'mean');
volume_aver_bin=groupsummary(m_big_up{:,3},depth_idx,'mean');
total_volume_bin=groupsummary(m_big_up{:,4},depth_idx,'mean');
ctd_rest_perim=m_big_up{:,7:26};
ctd_rest_bin=groupsummary(ctd_rest_perim,depth_idx,'mean');
ctd_rest_bin(:,16)=[];
date_bin={};
time_bin={};

for i = depth_group_uni'
        id_bin=find(depth_idx==i);
        date_bin=[date_bin;m_big_up(id_bin(1),5)];
        time_bin=[time_bin;m_big_up(id_bin(1),6)]; 
end
    
m_bin_ori_up=[table(depth_bin_uni,particle_no_bin,volume_aver_bin,total_volume_bin,table2cell(date_bin),table2cell(time_bin)),array2table(ctd_rest_bin)];
m_bin_ori_up=flip(m_bin_ori_up);
name_col={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time','Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2','days','hours','min','sec'};
m_bin_ori_up=renamevars(m_bin_ori_up,1:width(m_bin_ori_up),name_col);

%%
% find missing data in the table
% ori_ctd=readtable(['/home/hhuan006/CTD2021fall/' cs{2} '.csv'],1,0);

miss_dep=[];
dp_ctd=ori_ctd(:,16);
idx_max_dp_ctd=find(dp_ctd==max(dp_ctd));
dp_ctd_up=dp_ctd(idx_max_dp_ctd+1:length(dp_ctd));
ctd_up=ori_ctd(idx_max_dp_ctd+1:length(dp_ctd),:);
m_bin_up=m_bin_ori_up;
dp=table2array(m_bin_up(:,1));
name_col={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time','Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2','days','hours','min','sec'};
name_col_1={'Depth(m)','Particle no.','Average particle volume','Total particle volume','Date','Time'};
name_col_2={'Conductivity','Conductivity2','temperature','temperature2','pressure','oxy', 'oxy2', 'fluorescence[ug/l]','fluorescence[mg/m^3]','BeamAttenuation','BeamTransmission', 'PSU','PSU2','density', 'density2'};
name_col_3={'days','hours','min','sec'};
if dp(1)~=max(dp_ctd)-1
   for i = 1:max(dp_ctd)-dp(1)-1 
       id_ms_ctd=max(dp_ctd)-dp(1)-i;
       t=renamevars(array2table(ctd_up(id_ms_ctd,17:20)),1:4,name_col_3);
       column_ms=[array2table([ctd_up(id_ms_ctd,16),NaN,NaN,NaN,NaN,NaN],'VariableNames',name_col_1),array2table(ctd_up(id_ms_ctd,1:15),'VariableNames',name_col_2),t];

       m_bin_up=[column_ms;m_bin_up];
   end
end
if dp(length(dp))-dp(1)+1~=length(dp)
    for i = 1:length(dp_ctd_up)-1
        dp=m_bin_up{:,1};
        if dp(i)~=dp(i+1)+1
            
            miss_dep=[miss_dep;dp(i)-1];
            
            point_ms=dp(i)-1;
            id_ms_ctd=find(dp_ctd_up==point_ms);
            id_ms_img=min(find(dp==dp(i)));
            pre=m_bin_up(1:id_ms_img,:);
            post=m_bin_up(id_ms_img+1:height(m_bin_up),:);
            t=renamevars(array2table(ctd_up(id_ms_ctd,17:20)),1:4,name_col_3);
            column_ms=[array2table([point_ms,NaN,NaN,NaN,NaN,NaN],'VariableNames',name_col_1),array2table(ctd_up(id_ms_ctd,1:15),'VariableNames',name_col_2),t];

            m_bin_up=[pre;column_ms;post];
        end
    end
end

m_bin=[m_bin_down;m_bin_up];
end
save_path='/home/hhuan006/sum_2021fall/raw/smallthan1mm/T=29/';
% save_path='/home/hhuan006/sum_2021fall/all/T=29/sub100/';

if ~exist(save_path)
  mkdir(save_path)
end
writetable(m_bin,[save_path '/sum_' cs{2} '_' a{1} '.csv']);
toc



        

        
        
        

        