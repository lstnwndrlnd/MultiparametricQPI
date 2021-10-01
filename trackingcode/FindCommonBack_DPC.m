% This code 'FindCommonBack_DPC.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

clear all, close all
tic
wavelength = 624; %nm
% froot should point to raw phase data
froot = '';
% fstart should be the common beginning of all filenames in the experiment
fstart = '';
frames_used = 1000;

filelist = dir([froot, fstart, '*.mat']);
fileNames = char(sort_nat({filelist.name}'));
numf = length(fileNames);

%randomize file list
RNums = randperm(numf);
fileNames = fileNames(RNums,:);

B3D   = zeros(1200,1920,frames_used);

parfor frame = 1:frames_used
    fname = strtrim([froot, fileNames(frame,:)]);
    [D1,L1] = LoadSegment_short(fname, wavelength);


    BWfinal = L1>0;
    B3D(:,:,frame) = D1.*(~BWfinal);
end

dim = 3;
% the dimension along which mean is calculated

Btotal = mean(B3D,dim,'omitnan');
Ntotal = sum(~isnan(B3D),dim);
save([froot, 'Btotal'], 'Btotal', 'Ntotal');

toc