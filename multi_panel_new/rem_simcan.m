function [fill_edge,binary_close]=rem_simcan(diff,T)
%% This is the edge detection function constructed by Huanqing Huang (ODU)
% firstly, binarization image is obtained using global thresholding
% secondly, a canny edge detection is added to improve the particles

% use < when the background white and object black.
% use > when the object is white and background white
I=diff>T;

% bridge the gaps: sets 0-valued pixels to 1 if they have two nonzero neighbors that are not connected. 
b1=bwmorph(I, 'bridge');
% fill in the gap between pixels (4 connected by default)
b2=imfill(b1,'holes');

% filter out the particles with area < 3, it turned out 3 might be better
b3 = bwareaopen(b2, 3);
% figure
% imshow(binary1)
% create a structure element to do the morphological filters
se = strel('disk',1);
% use a dilation then erosion. This step is to remove the fluffy boundary
% of particles
binary_close = imclose(b3,se);
% b_ap=imcrop(binary_close,[192,458,10,10]);
% imwrite(b_ap,['D:\DQ\DQ\MS\research\pub\threshold algorithm\particle\glb_m_ap.tiff'])
% 
% figure
% imshow(binary_close)
% bw = imcrop(binary_close,[728,723,10,10]);
% figure
% imshow(bw)
% imshowpair(binary1, binary_close,'montage') % after opening and closing filters, there are more white pixels
% [~,edge_final]=Canny_no_sub(binary_close);
[~,edge_final]=Canny_detection(binary_close);

edge_final=bwmorph(edge_final, 'bridge');
%figure
% imshow(edge_final)
% fill in the pixels inside the edges
fill_edge=imfill(edge_final,'holes');
% %figure
% imshow(fill_edge)
% dilate then erode by the same kernel
% J= imclose(fill_edge, se);
% remove the edges in Canny detection
fill_edge(edge_final~=0)=0;
fill_edge=imfill(fill_edge,'holes');
% figure;imshow(fill_edge)
% ed_ap=imcrop(fill_edge,[192,458,10,10]);
% imwrite(ed_ap,['D:\DQ\DQ\MS\research\pub\threshold algorithm\particle\cann_ap.tiff'])

% filter out the particles with area < 3
fill_edge = bwareaopen(fill_edge, 3);
% fill again in case there are edges within the particles
fill_edge=logical(fill_edge);
% ca_ap = imcrop(fill_edge,[192,458,10,10]);
% figure;imshow(ca_ap)
fill_edge= imfill(binary_close | fill_edge,'holes');
% remove noises on the particle edge 
% fill_edge=bwmorph(fill_edge,'spur');



