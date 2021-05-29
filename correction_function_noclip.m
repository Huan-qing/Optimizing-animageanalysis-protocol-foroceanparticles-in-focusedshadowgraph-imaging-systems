function [corImage]=correction_function(blankImage, originalImage)

if size(blankImage,3)==3
    blankImage=rgb2gray(blankImage);
end
if size(originalImage,3)==3
    originalImage=rgb2gray(originalImage);
end

if ~isa(originalImage,'double')
    originalImage=im2double(originalImage);
end
% standard_region=imcrop(blankImage,[550, 450,99,99]);
% brightness = mean2(standard_region);

brightness = 207;
if brightness>1
    Gs=brightness/255;
else
    Gs=brightness;
end
% Convert the image format
% I = im2double(blankImage);
% Ig = im2double(originalImage);
%------ for deep sea images in 2019-2020, use the code below instead.-----
I = im2double(blankImage);
Ig = im2double(originalImage);

% Calculate the correction factor for each pixel point from blank image.
%-----------------------------------------------------------------------
Gain=Gs./I;
Gain(Gain==inf)=0;
% Clip any huge correction factor, because that will cause terrible
% noises!!!
% Gain(Gain>2)=2;
% Obtain the corrected blank image.
cor_blankimage=Gain.*I;
figure
imshow(cor_blankimage)
% Obtain the corrected gradient image.
cor_gradient=Gain.*Ig;
cor_gradient(cor_gradient==NaN)=0;
figure
imshow(cor_gradient)
corImage=imabsdiff(cor_blankimage,cor_gradient);

