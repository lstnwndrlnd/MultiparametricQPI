function [ D2, L, B1] = LoadSegment_Btotal_rollingball_BT474( fname, wavelength, Btotal)
% function to load and segment the data stored in fname
% phase data should be stored in fname.

% This code 'LoadSegment_Btotal_rollingball_BT474.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

% Btotal version loads Btotal to subtract it from the loaded image
%   This is usually used when running the tracking code

% TM 10/12/2020 Edited the load line to only load phase and use it directly
% This to save memory needed by half. Combined the steps to make D into 1
% step

load(fname,'Phase');

if isa(Phase,'single')
    D = double(Phase).*wavelength;
elseif isa(Phase,'int16')
    D = double(Phase)*(2*pi).*wavelength/65536;
else
    error('Error. \nInput must be a int16 or single, not a %s.',class(Phase))
end

[B1,M] = imagebackground_polyn_reduced2_v2(D-Btotal,10,10,10);    %compute image background
D1 = D-B1-Btotal;           %subtract background from data files
D2 = imtophat(D1, strel('ball', 100,100)); %40-50 worked well for MDA-MB-231
D2 = D2- mean(D2(M)); %make the background have 0 mean

DF = imfilter(D1, fspecial('gaussian', [12 12], 3));
L = imagesegment_aggressive_v4(DF); %segment image (detect distinct cell regions and disconnect connected regions)


end
