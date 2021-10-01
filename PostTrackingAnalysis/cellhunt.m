function [cells, cells_loc, cells_drug, cells_conc,cells_ID] = cellhunt(froot, minp, StartTime, EndTime)
%%% This code 'cellhunt.m' is copyright 2021
% Authors: Tarek E. Mustafa, Edward R. Polanco, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

minmass = 110;
%First load tracks and trackheader 
load([froot,'data_allframes.mat'],'tracks','tracksColHeaders');

% Use trackheader to find location of drugID, concentration, time and
% frameID columns locations
drugloc   = matches(tracksColHeaders,'condition_drugID');
concloc   = matches(tracksColHeaders,'condition_concentration (um)');
timeloc   = matches(tracksColHeaders,'Time (h)');
frameloc  = matches(tracksColHeaders,'Frame ID');
cellsloc  = matches(tracksColHeaders,'id');

% find the frame id that exist between starttime and endtime
target      = tracks(:,timeloc)>= StartTime & tracks(:,timeloc)<= EndTime;
frames      = unique(tracks(target,frameloc));
StartFrame  = min(frames);
EndFrame    = max(frames);

clearvars target frames


% find the cells id that exist between start frame and end frame
target      = tracks(:,frameloc)>= StartFrame & tracks(:,frameloc)<= EndFrame;

% the next line will find cells ID between frames then
cellsID       = unique(tracks(target,cellsloc));
Ncount        = histc(tracks(target,cellsloc),cellsID); 
targetcells   = cellsID(Ncount>= minp);

newtracks  = tracks(target,:);


remove= zeros(1,length(targetcells)); % Empty array to track how many cells were removed
mass  =  newtracks(:,3);
IDs   = newtracks(:,5);

% can be parallelized with parfor
parfor i = 1:length(targetcells)
    y = IDs==targetcells(i);
    m = mean(mass(y));
    if m <= minmass
        remove(i) = 1;
    end
end
z   = find(remove);
for n = 1:length(z)
    targetcells(z(end-n+1)) = [];
end
    


cells         = cell(1,length(targetcells));
cells_loc     = zeros(1,length(targetcells));
cells_drug    = zeros(1,length(targetcells));
cells_conc    = zeros(1,length(targetcells));
cells_ID      = zeros(1,length(targetcells));

% can be parallelized with parfor
parfor z = 1:length(targetcells)
    
    target_tracks  = newtracks(newtracks(:,5) ==targetcells(z),:);
%     cells{z}      =[target_tracks(:,timeloc) target_tracks(:,3)];
    cells{z}      =[target_tracks(:,timeloc) medfilt1(target_tracks(:,3),5,'truncate')];
    cells_loc(z)  = target_tracks(1,8);
    cells_drug(z) = target_tracks(1,drugloc);
    cells_conc(z) = target_tracks(1,concloc);    
    cells_ID(z)   = target_tracks(1,5); 
end

end