function [ DF, L ] = LoadSegment_Btotal( fname, wavelength, Btotal )
% function to load and segment the data stored in fname
% phase data should be stored in fname.Phase

% This code 'LoadSegment_Btotal.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.


Loaded = load(fname);
D = double(Loaded.Phase);
D = D.*wavelength;


B = imagebackground_polyn_reduced2(D,10,50,10);    %compute image background
D = ((D-B-Btotal));           %subtract background from data files
DF = medfilt2(D); % median filter background corrected phase image
L = imagesegment_aggressive(DF); %segment image (detect distinct cell regions and disconnect connected regions)


end

