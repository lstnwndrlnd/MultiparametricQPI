function [M] = backgroundmask(I)
% function [M] = backgroundmask(I)
% function to mask the background of an image, I
% input: I, the grayscale image to find the background of
% output: B, the background of I
% method: find 'objects' in I, mask them from the image, fit the remaining
% pixels with an ellipse, subtract that from the image

% This code 'backgroundmask.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.


[junk threshold] = edge(I, 'sobel');
fudgeFactor = .4;
%fudgeFactor = 0.8;
BWs = edge(I,'sobel', threshold * fudgeFactor);

%step 2, dilate the image
se90 = strel('line', 4, 90);
se0 = strel('line', 4, 0);
BWsdil = imdilate(BWs, [se90 se0]);
BWdfill = imfill(BWsdil, 'holes');

BWf = imerode(BWdfill, strel('disk', 2));

M = bwareaopen(BWf,50); %remove all areas of fewer than 50 pixels


end



