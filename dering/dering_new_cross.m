% close all
clear
path='/home/hhuan006/dering/';
% path='D:\DQ\DQ\MS\research\optical experiment\blur\ring_exp\apend\';
FileList   = dir(fullfile(path, '*.tif*'));
tic
for i= 12%:length(FileList) %  good 14 12
    clear edge_coor
    image_name = FileList(i).name;
    grayImage=im2double((imread(strcat(path,image_name))));
dering_img=grayImage;

obj_binary=rem_simcan(grayImage, 15/255);
main_edge=bwmorph(obj_binary,'remove');
rep_bright=median(grayImage((grayImage>graythresh(grayImage))~=0));
% figure;
% imshow(obj_binary)
 
 % find the cross-section
 num_edge=sum(main_edge(:));
 idx_edge=find(main_edge==1);
 [edge_coor(:,2),edge_coor(:,1)]=find(main_edge~=0); % the x and y coordinate are reversed here
 figure
 imshow( grayImage)
 
    % Method 2 : to find the major axis and minor axis, then draw perpendicular line
    
dist_min=inf;
dist_max=-inf;

for p = 1: num_edge

    coor_edge=[edge_coor(p,:)];
    for q =1: num_edge
        coor_sk=[edge_coor(q,:)];
        if ~isequal(coor_edge,coor_sk)
        dist_one=norm(coor_edge-coor_sk);
        
        if dist_one>dist_max
            dist_max=dist_one;
            major_ax(1,:)=coor_edge;
            major_ax(2,:)=coor_sk;
        end
        
        end
    end
end

obj_dist_max=0;
% get the extended line of major axis and its coordinate
[major_cross]=drawLine_perpendicular(size(grayImage),major_ax);
line([major_cross(1,1) major_cross(2,1)],[major_cross(1,2) major_cross(2,2)],'color','yellow')
[cx_major,cy_major,~,~,~]=improfile(grayImage,[major_cross(1,1) major_cross(2,1)],[major_cross(1,2) major_cross(2,2)]);

% find minor axis for each coordinate on major axis (slope not change)
slope_major=(major_cross(1,2)-major_cross(2,2))/(major_cross(1,1)-major_cross(2,1));
midpoint=[mean(major_ax(:,1)),mean(major_ax(:,2))];
slope_minor=-1/slope_major;

for p = 1: length(cx_major)

    [minor_cross]=drawLine_perpendicular(size(grayImage),[cx_major(p),cy_major(p)],slope_minor);
    line([minor_cross(1,1) minor_cross(2,1)],[minor_cross(1,2) minor_cross(2,2)],'color','red')

% figure;improfile(grayImage,[x1 x2],[y1 y2])
[cx,cy,c,xi2,yi2]=improfile(grayImage,[minor_cross(1,1) minor_cross(2,1)],[minor_cross(1,2) minor_cross(2,2)]);
coor_line=[cx,cy];

ind_max=find(c==max(c));

if max(c)<0.1 || ind_max(1)==1 || ind_max(end)==length(c)  % filter out the signal that are not meaningful
    continue
else
[sig_max,loc_max]=findpeaks(c);
[sig_min,loc_min]=findpeaks(-c);

pk_num=0;
ring_coor=[0,0];
local_rep=max(c);
ind_range=find(c>=local_rep*0.6);
% decide the obj signal range
if sum(c(ind_range)>local_rep*0.1)~=length(ind_range) % if there is a gap in the signal range
    ind_min=find(c(ind_range)<=local_rep*0.1);
    if length(ind_min)==1
        if ind_min>floor(length(ind_range)/2)
            ind_range=ind_range(1:ind_min);
        else 
            ind_range=ind_range(ind_min:end);
        end
    else
        series_gap=[1,ind_min,length(ind_range)];
        dist_gap=diff(series_gap);
        ind_main=find(dist_gap==max(dist_gap));
        ind_range=ind_range(ind_main-1:ind_main);
    end
end
% ind_range=ind_range(diff(ind_range)==1);
% ind_big=loc_max(find(sig_max>=local_rep));
sr=[1:length(loc_max)];
ind_globmax=sr((ismember(loc_max,ind_range))); % find peaks that is greater than 
if ~isempty(ind_globmax)
    max_sig=max(sig_max);
    obj_left=loc_min(loc_min<min(ind_range));
     
    if isempty(obj_left)
        obj_1=1;
    else
        obj_1=obj_left(end);
    end
    obj_right=loc_min(loc_min>max(ind_range));
    
    if isempty(obj_right)
        obj_2=length(c);
    else
        obj_2=obj_right(1);
    end
    
    obj_dist_coor=[obj_1(1),obj_2(end)];
    obj_dist=norm(coor_line(obj_dist_coor(1),:)-coor_line(obj_dist_coor(2),:));
    if obj_dist>obj_dist_max
        obj_dist_max=obj_dist;
        obj_dist_max_coor=obj_dist_coor;
    end
        
    for i_pk=1:floor(length(sig_max)/2)-1
      if ind_globmax(1)-i_pk>0
       
        if sig_max(ind_globmax(1)-i_pk)>0.1*max_sig && all(diff(sig_max(ind_globmax(1)-i_pk:ind_globmax(1)))>0) % signal stronger than 10 % and peak monotonic 
           if ind_globmax(1)-i_pk>1
            loc_r1=loc_max(ind_globmax(1)-1-i_pk)+find(c(loc_max(ind_globmax(1)-1-i_pk):loc_max(ind_globmax(1)-i_pk))==min(c(loc_max(ind_globmax(1)-1-i_pk):loc_max(ind_globmax(1)-i_pk))))-1;
            ring_coor(1,1)=loc_r1(1);
           else
            ring_coor(1,1)=1;
           end
        end
      end
      if ind_globmax(end)+i_pk<=length(sig_max)
          if sig_max(ind_globmax(end)+i_pk)>0.1*max_sig && all(diff(sig_max(ind_globmax(end):ind_globmax(end)+i_pk))<0) % signal stronger than 10 % and monotonic
            if ind_globmax(end)+i_pk<length(sig_max)
              loc_r2=loc_max(ind_globmax(end)+i_pk)+find(c(loc_max(ind_globmax(end)+i_pk):loc_max(ind_globmax(end)+i_pk+1))==min(c(loc_max(ind_globmax(end)+i_pk):loc_max(ind_globmax(end)+i_pk+1))))-1;
              ring_coor(1,2)=loc_r2(end);
            else 
                ring_coor(1,2)=length(c(~isnan(c)));
            end
          end
      end
      end
    
 end

if ring_coor(1)==0
    
    ring_ratio(p,1)=0;
    
else
 
    ring_dist1=norm(coor_line(obj_dist_coor(1),:)-coor_line(ring_coor(1),:));
    ring_ratio(p,1)=ring_dist1/obj_dist;
    
end

if ring_coor(2)==0
    
    ring_ratio(p,2)=0;
    
else
 
    ring_dist2=norm(coor_line(obj_dist_coor(2),:)-coor_line(ring_coor(2),:));
    ring_ratio(p,2)=ring_dist2/obj_dist;
    
end
end

% dering process
for idd = 1:2
    if ring_ratio(p,idd) ~= 0
        pt1 = ring_coor(idd);
        pt2 = obj_dist_coor(idd);
        v_1(idd)=c(pt1);
        v_2(idd)=c(pt2);
    end
end 
% dering to the artifact in the left

% dering_img(1:obj_dist_coor(2),1:obj_dist_coor(1))=0;
for  f=1:length(coor_line)
    
if f<obj_dist_coor(1)||f>obj_dist_coor(2) % if the coordinate index in the range out of obj edge
    coor_one=round(coor_line(f,:));
    coor_one(coor_one==0)=1;
    dering_img(coor_one(2),coor_one(1))=0;
    f=f+1;
end
end
% figure;imshow(dering_img)
% f=1;
% dering to the right
 end
            
mean(ring_ratio(:))  
          
% im_dering=grayImage;
% im_dering(newb==1)=0;
ori_binary=rem_simcan(grayImage,20/255);
gray_dering=imgaussfilt(dering_img,0.7);
dering_binary=rem_simcan(gray_dering,20/255);
% figure;imshow(grayImage)
% 
% figure;imshow(ori_binary)
% 
% figure;imshow(gray_dering)
% 
% figure;imshow(dering_binary)

imwrite(grayImage,['/home/hhuan006/dering/example/ori' num2str(i) '.bmp'])
imwrite(ori_binary,['/home/hhuan006/dering/example/bi_ori' num2str(i) '.bmp'])
imwrite(gray_dering,['/home/hhuan006/dering/example/der' num2str(i) '.bmp'])
imwrite(dering_binary,['/home/hhuan006/dering/example/bi_der' num2str(i) '.bmp'])

% figure;imshow(dering_img>rep_bright/255)

% figure;imshow(grayImage)
% improfile(grayImage,[48,85],[11,66])
% [cx,cy,c,xi2,yi2]=improfile(grayImage,[90 95],[1 100]);
% improfile(grayImage,[90 95],[1 100])
% % figure,imshow(B)
% figure;imshow(mask_gt50)
% figure;coor_line
%  imshow(main_edge)
%  figure
%  imshow(labeloverlay(grayImage,skeout,'Transparency',0))
% figure;
% imshow(obj_edge)
% figure;imshow(grow_fill)
% figure;imshow(newb)

end
toc