% function [ D, L, B ] = LoadSegment_Btotal_short_rollingball( fname, wavelength, Btotal )
function [ D2, L, B1] = LoadSegment_Btotal_rollingball_MDA( fname, wavelength, Btotal )
% function to load and segment the data stored in fname
% phase data should be stored in fname.Phase

% This code 'LoadSegment_Btotal_rollingball_MDA.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

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
D = D-B1-Btotal;           %subtract background from data files
D2 = imtophat(D, strel('ball', 70,70)); %40-70 worked well for MDA-MB-231
D2 = D2- mean(D2(M)); %make the background have 0 mean

DF1 = medfilt2(D2); % median filter background corrected phase image
DF2 = imfilter(DF1, fspecial('gaussian', [12 12], 3));
L = imagesegment_aggressive_v3(DF2); %segment image (detect distinct cell regions and disconnect connected regions)


end
