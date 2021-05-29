function [fill_edge]=rem_simcan(diff,thre)

% use < when the background white and object black.
% use > when the object is white and background white
binary=diff>thre;

%figure
% imshow(binary)
% bridge the gaps: sets 0-valued pixels to 1 if they have two nonzero neighbors that are not connected. 
binary=bwmorph(binary, 'bridge');
% fill in the gap between pixels (4 connected by default)
binary1=imfill(binary,'holes');
% filter out the particles with area < 3, it turned out 3 might be better
binary1 = bwareaopen(binary1, 3);
% create a structure element to do the morphological filters
se = strel('disk',1);
% use a dilation then erosion. This step is to remove the fluffy boundary
% of particles
binary_close = imclose(binary1,se);
%%figure
% imshow(binary_close)
% %%figure
% imshowpair(binary1, binary_close,'montage') % after opening and closing filters, there are more white pixels
% [~,edge_final]=Canny_no_sub(binary_close);
[~,edge_final]=Canny_no_sub(binary_close);
edge_final=bwmorph(edge_final, 'bridge');
%figure
% imshow(edge_final)
% fill in the pixels inside the edges
fill_edge=imfill(edge_final,'holes');
% %figure
% imshow(fill_edge)
% dilate then erode by the same kernel
J= imclose(fill_edge, se);
% remove the edges in Canny detection
fill_edge(edge_final~=0)=0;
fill_edge=imfill(fill_edge,'holes');
% filter out the particles with area < 3
fill_edge = bwareaopen(fill_edge, 3);
% fill again in case there are edges within the particles
fill_edge=logical(fill_edge);
%figure;imshow(fill_edge)
fill_edge= imfill(binary_close | fill_edge,'holes');
% remove noises on the particle edge 
% fill_edge=bwmorph(fill_edge,'spur');



