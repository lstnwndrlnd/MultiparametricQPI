% This code 'tracking_DPC.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

% script to run cell tracking codel-ko
% TAZ 5/18/10, modified 7/21/10 to automatically detect file names and to
% remove drift due to the entire frame moving
% first section: define image processing parameters (file locations,
% processing options, etc.)

% second section: load images and detect cells
% third section: track cells

clear all; close all;
tic
%%
froot = '';
fstart = '';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% comment or uncomment the appropriate lines to choose cell line and drug
% panel

% panel = 1 if drug panel A
panel = 2;
% panel = 2 if drug panel B
% panel = 2;

% cellLine = 1 if MCF7
cellLine = 1;
% cellLine = 2 if MDA-MB-231
% cellLine = 2;
% cellLine = 3 if BT-474
% cellLine = 3;

wavelength = 624; %nm
pxlsize = 5e-4; %mm/pixel, for 10x
savefolder = froot;
% savefolder = 'F:\Eddie\DrugResponseExperiments\MDA-MB-231\PanelB\ExploratoryExperiment_090120\datasubsets\fulldata\';

%%% First, define which files the script will work on
fext = '.mat'; %file extension
overwrite = 0; %set to 1 to enable overwrite of pre-stored data files

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Second, define image analysis parameters

%%% define min and max area and mean intensity (MI) thresholds
%%% only objects which fall between these values will be counted as "cells"
%%% These parameters should be adjusted for each sample to capture the
%%% objects of interest. See Figure 11 to evaluate where these values fall
%%% relative to the properties of the image. See figures 12 and 13 to see
%%% which objects in the first and last frames are counted as "cells"
% minAreathresh = 500;
% 10um dataset

% For MDA use 1000 30000 100 800
% For MCF7  1000 222222 50 800
% MFC7
if (cellLine == 1 || cellLine == 3)
    minAreathresh = 1000; 
    maxAreathresh = 222222; 
    minMIthresh = 50;
    maxMIthresh = 800;
elseif (cellLine == 2)
    % MDA
    minAreathresh = 900; 
    maxAreathresh = 100000; 
    minMIthresh = 80;
    maxMIthresh = 800;
end
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Third, define tracking software parameters %%%%%%%%%%%%%%%%%

minpathlength = 20; %min path length to use in plotting results. only paths
%                    of this length or longer will be displayed. this does
%                    not affect the tracking software (tracks shorter than
%                    minpathlength will still be computed and stored)

massfact = 0.25;   %factor to multiply mass by in tracking step. use this to
%account for differences in how much the cell moves vs. how
%much mass changes over time

%%% tracking parameters below affect the tracking software itself
max_disp = 40; %max displacement for particle tracking
%               max_disp is an estimate of the maximum distance that a
%               particle would move in a single time interval. It should be
%               set to a value somewhat less than the mean spacing between
%               the particles
param.mem = 1;  %this is the number of time steps that a particle can be
%               'lost' and then recovered again.  If the particle reappears
%               after this number of frames has elapsed, it will be
%               tracked as a new particle. The default setting is zero.
%               this is useful if particles occasionally 'drop out' of
%               the data.
param.dim = 3; %number of dimensions of coordinate data to track. If you
%               set this value to 2, then it will track based on position
%               in x and y. If param.dim = 3, then the software will track
%               based on mass as well.
param.good = 0; %set this keyword to eliminate all trajectories with
%                fewer than param.good valid positions.  This is useful
%                for eliminating very short, mostly 'lost' trajectories
%                due to blinking 'noise' particles in the data stream.
param.quiet = 1; %set this keyword to 1 if you don't want any text
%                 displayed while the tracking algorithm is running


%%

%%
[LocList, numLoc] = getloclist(froot, fstart, fext);
%pre-processing and variable initialization before loop begins:
Loc = 1;
filelist = dir([froot, fstart, char(LocList(Loc)), '_*', fext]);
fileNames = char(sort_nat({filelist.name}'));
numf = length(fileNames);

fname = strtrim([froot, fileNames(1,:)]);

load([froot, 'Btotal'], 'Btotal');
if (cellLine == 1)
    [D1,~] = LoadSegment_Btotal_rollingball_MCF7(fname,wavelength, Btotal);
elseif (cellLine == 2)
    [D1,~] = LoadSegment_Btotal_rollingball_MDA(fname,wavelength, Btotal);
elseif (cellLine == 3)
    [D1,~] = LoadSegment_Btotal_rollingball_BT474(fname,wavelength, Btotal);
end

time0 = LoadTime(fname);

D1s = zeros([size(D1),numLoc], 'single');

%%
%loop over all locations
parfor Loc = 1:864 %parfor                                                                                                                                                                                  
    if ~exist([savefolder, 'Tdata', num2str(Loc), '.mat'],'file') || overwrite
%         Loc = 1;
        filelist = dir([froot, fstart, char(LocList(Loc)), '_*', fext]);
        fileNames = char(sort_nat({filelist.name}'));
        numf = length(fileNames(:,1)); % track all frames at position

        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% grab first frame for analysis and detection of the correct cell
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        fname = strtrim([froot, fileNames(1,:)]);

        if (cellLine == 1)
            [D1,L1] = LoadSegment_Btotal_rollingball_MCF7(fname,wavelength, Btotal);
        elseif (cellLine == 2)
            [D1,L1] = LoadSegment_Btotal_rollingball_MDA(fname,wavelength, Btotal);
        elseif (cellLine == 3)
            [D1,L1] = LoadSegment_Btotal_rollingball_BT474(fname,wavelength, Btotal);
        end
        
        %preallocate variables for speed
        yshift_store = zeros(numf);
        xshift_store = zeros(numf);
        t_stored = zeros(numf);
        
        % these can be commented out to save memory, if you do that then
        % make sure to comment out the lines that use them belows
        D_stored = zeros([size(D1),numf], 'int16');
        L_stored = zeros([size(D1),numf], 'uint16');
        
        
        D1s(:,:,Loc) = single(D1);
        
        %%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% loop through first numf file names stored in fnum and store analysis
        %%% results
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tt = 1; %initialize tt, the index of the tracking array
        T_array = [];
        
        xshift_old = 0;
        yshift_old = 0;
        D_old = D1;
        segID = [];
        labelID = [];
        
        for jj = 1:numf
            
            fname = strtrim([froot, fileNames(jj,:)]);
            disp(fname)

            if (cellLine == 1)
                [D,L] = LoadSegment_Btotal_rollingball_MCF7(fname,wavelength, Btotal);
            elseif (cellLine == 2)
                [D,L] = LoadSegment_Btotal_rollingball_MDA(fname,wavelength, Btotal);
            elseif (cellLine == 3)
                [D,L] = LoadSegment_Btotal_rollingball_BT474(fname,wavelength, Btotal);
            end

            timen = LoadTime(fname);
            time = (datenum(timen)-datenum(time0)).*24; %store time in hours
            
            [V, M, A, MI, P, SF] = imageprops_SF(L, D, pxlsize); %compute image properties based on the regions stored in L
            % outputs of imageprops are in same order as labels in label
            % matrix, L. Therefore, we just need labelID to be 1:max(L(:))
            
            if std(D(:))~=0 %skip if blank image
                [yshift, xshift] = CorrShift(D_old,D); %find average shift between current frame and first frame
                yshift = yshift+yshift_old;
                xshift = xshift+xshift_old;
                D_old = D;
                xshift_old = xshift;
                yshift_old = yshift;
                
                yshift_store(jj,Loc) = yshift;
                xshift_store(jj,Loc) = xshift;
                
                D_stored(:,:,jj) = int16(D(:,:)/2/pi*65536/wavelength);
                L_stored(:,:,jj) = uint16(L(:,:));
                t_stored(jj,Loc) = time;
                
                %next, loop through all items identified in V and find only the ones
                %which meet area and mean intensity requirements
                for ii = 1:length(V)
                    %first, check that 1) there is something at index ii, 2) that
                    if(~isnan(P(ii).Centroid(1)) && A(ii) > minAreathresh && A(ii) < maxAreathresh && MI(ii) > minMIthresh && MI(ii) < maxMIthresh)
                        T_array(tt,1:2) = P(ii).Centroid; %store position in first two columns of T_array
                        T_array(tt,1:2) = T_array(tt,1:2) - [xshift, yshift]; %remove shift due to movement of the entire frame
                        T_array(tt,3)   = M(ii);          %store mass in third column
                        T_array(tt,4)   = time;           %store time from first frame in seconds in 4th column
                        T_array(tt,5)   = A(ii); %store area in fifth column
                        T_array(tt,6)   = SF(ii); %store shape factor in sixth column
                        labelID(tt) = ii;
                        tt = tt+1;                        %increment T_array index
                    end
                end
            end
        end
        segID = [segID; labelID];
        
        if ~isempty(T_array) && sum(T_array(:,4) ~= T_array(1,4))>0
            %%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Cell tracking using the track function
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % T_array starts as [x, y, m, t, A, SF]
            T_array(:,3) = T_array(:,3).*massfact; %change mass weighting in T_array
            T_array = sortrows(T_array, 4); %sort T_array based on time vectors
            minTx =  min(T_array(:,1));
            T_array(:,1) = T_array(:,1) -minTx +1; %make sure all x positions are positive\
            minTy =  min(T_array(:,2));
            T_array(:,2) = T_array(:,2) -minTy +1; %make sure all y positions are positive, new with rev6
            % move time to last column, T_array will now be [x, y, m, A, SF, t]
            T_array_temp = [T_array(:,1:3), T_array(:,5:6), segID', T_array(:,4)];
            
            tracks = track(T_array_temp,max_disp,param);
            tracks0 = tracks;
            
            tracks_temp = [tracks(:,1:3), tracks(:,7:8), tracks(:,4:5)];
            sortedSeg = tracks(:,6);
            tracks = tracks_temp;
            
            T_array(:,3) = T_array(:,3)./massfact; %adjust mass weighting back to the way it was
            T_array(:,1) = T_array(:,1) + minTx -1; %set all x positions back to the way they were
            T_array(:,1) = T_array(:,2) + minTy -1; %set all y positions back to the way they were
            tracks(:,3) = tracks(:,3)./massfact; %adjust for mass weighting in the tracks array as well
            tracks(:,1) = tracks(:,1) +minTx -1; %set all x positions back to the way they were
            tracks(:,2) = tracks(:,2) +minTy -1; %set all y positions back to the way they were
        else
            % initialize tracks and sortedSeg in case condition not met
            % above
            tracks = [];
            sortedSeg = [];
        end
        
        % for each location, insert each location number
        tracks(:,8) = Loc;
        times = tracks(:,4);
        
        uniqueTimes = t_stored(:,Loc);
        
        % numf is also the length of uniqueTimes, so for each unique time:
        for ii = 1:numf
            % associate frame number with timestamp on the image
            currentTime = tracks(:,4) == uniqueTimes(ii);
            % insert frame number
            tracks(currentTime,9) = ii;
        end
        
        if(panel == 1)
            [drugID, conc, ~, ~, ~] = getCondition_PanelA_v1(Loc);
        elseif(panel == 2)
            [drugID, conc, ~, ~, ~] = getCondition_PanelB_v1(Loc);
        end
            
        
        sizeTracks = size(tracks);
        numRows = sizeTracks(1);
        % check to make sure that there is something in segID 
        if isempty(segID) || numRows == 1
            % if empty there should only be one row in tracks, set 10
            % element to zero
            tracks(:,10) = 0;
        else
            % if it isn't empty then the size of these arrays should match5
            tracks(:,10) = sortedSeg;
        end
        tracks(:,11) = drugID;
        tracks(:,12) = conc;
        
        
        % save D_stored and L_stored in a separate .mat file to save memory
        
        % save(['data', row, col, '.mat'], 'D_stored', 'L_stored') %
        % line above is an alternative to the two lines below (I think the
        % parsave() is for saving in parallel to speed up processing.
        parsave([savefolder 'data', num2str(Loc), '.mat'], D_stored, 'D_stored', 0)
        parsave([savefolder 'data', num2str(Loc), '.mat'], L_stored, 'L_stored', 1)
        
        %save tracks data and clear vector so that code can be parallelized
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], tracks, 'tracks', 0)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], T_array, 'T_array', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], xshift_store, 'xshift_store', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], yshift_store, 'yshift_store', 1)
        parsave([savefolder 'Tdata', num2str(Loc), '.mat'], t_stored, 't_stored', 1)
        
    end %close if statement looking for stored data
end %close rows for loop

clear tracks T_array xshift_store yshift_store t_stored D_stored L_stored T_array_temp tracks_temp
toc
disp('End of parfor loop')
tic
datetime('now')
%%
tic
%reconstitute tracks vector
tracks_comp = [];
xshift_store_c = [];
yshift_store_c = [];
t_stored_c = [];
ii_stored = [];
maxindex = 0;
Loc_stored_c = [];
% for Loc = 1:numLoc
for Loc = 1:864
    load([savefolder, 'Tdata', num2str(Loc), '.mat']);
    if ~isempty(tracks)
        tracks(:,5) = tracks(:,5)+maxindex;
    end
    
    tracks_comp = [tracks_comp; tracks];
    xshift_store_c = [xshift_store_c; xshift_store(:,Loc)];
    yshift_store_c = [yshift_store_c; yshift_store(:,Loc)];
    t_stored_c = [t_stored_c; t_stored(:,Loc)];
    ii_stored = [ii_stored, 1:length(t_stored(:,Loc))'];
    Loc_stored_c = [Loc_stored_c, (1:length(t_stored(:,Loc))').*0 + Loc];
    
    if ~isempty(tracks_comp)
        maxindex = max(tracks_comp(:,5));
    end
    
    clear tracks xshift_store yshift_store t_stored
    
end

Loc_stored = Loc_stored_c;
tracks = tracks_comp;
xshift_store = xshift_store_c;
yshift_store= yshift_store_c;
t_stored = t_stored_c;
ii_stored = ii_stored';
clear tracks_comp xshift_store_c yshift_store_c t_stored_c maxindex Loc_stored_c
if length(tracks) > 10
    T0=min(tracks(:,4)); %find time of first image in the set
    tracks(:,4) = (tracks(:,4)-T0);
    t_stored = t_stored - T0;
end
toc

%%
%%% added these two lines to remove superfluous images and save this data %%%
tic
clear D D1 D1s D_old L L1 Btotal
if(panel == 1)
    [drugID, conc, ~, ~, ~] = getCondition_PanelA_v1(Loc);
elseif(panel == 2)
    [drugID, conc, ~, ~, ~] = getCondition_PanelB_v1(Loc);
end 

tracksColHeaders = {'X', 'Y', 'Mass (pg)', 'Time (h)', 'id', 'Area', 'shape factor', ...
                            'Location ID', 'Frame ID', 'segmentLabel', 'condition_drugID', ...
                            'condition_concentration (um)'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~exist([savefolder 'data_allframes']) || overwrite
    save([savefolder 'data_allframes'])
end

toc
disp('data_allframes is created')
datetime('now')

