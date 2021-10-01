function [L,BWdfill] = imagesegment_aggressive_watershed(I)
% function to segment an image of cells using the matlab help file
% procedure followed by breaking up connected cells using a watershed image
% transform
% input: I, the grayscale image to segment
% output: L, the labeled, segmented image
% modified heavily for DPC image analysis

% This code 'imagesegment_aggressive_watershed.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.


BWs = edge(I,'sobel');

se0 = strel('line', 4, 0);
se90 = strel('line', 4, 90);
se = [se0 se90];
BWsdil1 = imdilate(BWs, se);

BWclear = imclearborder(BWsdil1); %remove regions smaller than 300 pixels

se0 = strel('line', 2, 0);
se90 = strel('line',2,90);
se = [se0 se90];
BWer1 = imerode(BWclear, se);
BWer2 = imerode(BWer1, se);

BWopen = bwareaopen(BWer2, 300);
% smooth image
seD = strel('diamond',3);

seCirc = strel('disk', 10);
BWdil2 = imdilate(BWopen,seCirc);
BWdil3 = imdilate(BWdil2,seCirc);
% BWclear = imclearborder(BWopen);
BWfinal = imfill(BWdil3,'holes');

Igr = mat2gray(I);
thr = 0.15; 
mask_em = imextendedmax(Igr, thr); 
mask_em = imfill(mask_em, 'holes');
%Next step: complement the image so that the peaks become valleys. We do
%this because we are about to apply the watershed transform, which
%identifies low points, not high points. 
I_c = imcomplement(I);
%Next: modify the image so that the background pixels and the extended
%maxima pixels are forced to be the only local minima in the image. 
I_mod = imimposemin(I_c, ~BWfinal | mask_em);
%now, compute watershed transform
L = watershed(I_mod);