function [L] = imagesegment_aggressive(I)
%function [L] = imagesegment_aggressive_v3(I)
%function to segment an image of cells using a watershed image
%transform
%input: I, the grayscale image to segment;
%output: L, the labeled, segmented image
%modified heavily for DPC image analysis

% This code 'imagesegment_aggressive.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

If1 = imfilter(I, fspecial('sobel'));
If2 = imfilter(I, fspecial('sobel')');
BWs = abs(If1)>40 | abs(If2)>40;


se90 = strel('line', 8, 90);
se0 = strel('line', 8, 0);

BWsdil = imdilate(BWs, [se90 se0]);

BWdfill = imfill(BWsdil, 'holes');

seD = strel('diamond',1);

BWfinal = imerode(BWdfill,seD);
BWfinal = bwareaopen(BWfinal, 1000); %remove regions smaller than 300 pixels

L1 = imclearborder(BWfinal);

L = bwlabel(L1);
