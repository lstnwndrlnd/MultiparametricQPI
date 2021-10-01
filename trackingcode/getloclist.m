function [LocList, numlocs] = getloclist(froot, fstart, fext)
% function to get the list of location strings for all files within a given
% directoy
% will skip any filenames starting in the charcter '%'

% This code 'getloclist.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

filelist = dir([froot, fstart, '*', fext]);
fileNames = char(sort_nat({filelist.name}'));
numf = length(fileNames);

strs = cell(numf-1,1);

for ii = 1:numf-1
    S = textscan(fileNames(ii,:), '%[QPM10x]%[051921]%[pos]%s', 'delimiter', '_');
    strs{ii} = char(S{4});
end

LocList = unique(strs);
LocList = LocList(~strcmp(LocList, ''));
LocList = sort_nat(LocList);

numlocs = length(LocList);
