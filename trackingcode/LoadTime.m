function [time] = LoadTime(fname)
% function [time] = LoadTime(fname)
% input: filepath including name of file with extension (matfile)
% output: reads the time stored in the mat file structure, time is returned
% in seconds
% alternate use: can be readily modified to read the timestamp off of the
% metadata of an image file (.jpg, .tiff, etc).

% function to load the time associated with a given filename
% this time will be a datestring, as returned by the dir function
% fname is the .mat file name associated with the .tif file of the original
% data, so time is the modification time of the original tif file.

% This code 'LoadTime.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license lets you to remix, adapt, and build 
% upon this work, as long as they credit the authors and distribute the
% modified code under the same license as the original.

% time is in seconds.
Loaded = load(fname);
% uncomment next line to read time of image metadata
% time = Loaded.TimeStamp.datenum;
% uncomment next line to extract the time variable from the .mat file
time = Loaded.timestamp;

end

