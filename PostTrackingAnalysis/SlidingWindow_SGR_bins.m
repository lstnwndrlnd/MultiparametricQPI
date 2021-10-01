%%% This code 'SlidingWindow_SGR_bins.m' is copyright 2021
% Authors: Tarek E. Mustafa, Edward R. Polanco, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.
%%
clear all, close all
%%
% This code takes length of the experiment and Window length to make
% multple plots.
% This code can be modified to use technical replicates (Same drugs & cell
% line experiments) instead of slices
% The code will always run from time = 0 to required time.

% froot should point to folder containing data_allframes.mat
froot = './';
% uncomment next line and last line if you would like to save results automatically
% savename = 'slidingwindowresults';

% ethanolcontrol variable points to 4-HT condition because that is the one that used 
% ethanol as a control code assumes ethanol is drug label 1

% uncomment the two lines below either `Panel A:` or `Panel B:` but not both
% Panel A:
% druglabels = {'EtOH','DMSO','4-HT','Lapt','Fulv','Doxo','5Fu'}; % For Panel A
% ethanolcontrol   = 3;

% Panel B:
druglabels = {'Untreated', 'DMSO', 'Palbociclib', 'Vinblastine', 'Carboplatin', 'Docetaxel'};
ethanol control = 0;

ExpFrames        = 216; % no of frames per experiment
SingleWindow     = 73;  % no of frames for each interval
filtersize       = 5;   % median filter element size
minmass          = 130; % minimum mean mass for cell to be considered a debris
minp             = 25;  % minimum path length to be considered a cell
binningRate      = 8;   % Number of intervals to bin together

tic
for i=1:ExpFrames - SingleWindow +1
    frame(1,i) = i;
    frame(2,i) = SingleWindow +i -1;
end

b_cell            = cell(1,i); 
cells_loc_cell    = cell(1,i);
cells_drug_cell   = cell(1,i);
cells_conc_cell   = cell(1,i);
cells_ID_cell     = cell(1,i);
b_raw_cell             = cell(1,i);
a_raw_cell             = cell(1,i);
cells_intmass_cell     = cell(1,i);
cells_meanmass_cell    = cell(1,i);
cells_finalmass_cell   = cell(1,i);

load([froot,'data_allframes.mat'],'tracks','tracksColHeaders');
% can be parallelized using parfor
parfor i = 1:ExpFrames - SingleWindow + 1
    [b_cell{i},cells_loc_cell{i},cells_drug_cell{i},cells_conc_cell{i},cells_ID_cell{i},...
        b_raw_cell{i},a_raw_cell{i},cells_intmass_cell{i},cells_meanmass_cell{i},cells_finalmass_cell{i}] = ...
        sgrmaker_mean(tracks,tracksColHeaders,frame(1,i),frame(2,i),minmass,minp);
end

% Convert cell arrays into long vectors
% First make the vector with first cell
b = b_cell{1};
b_raw = b_raw_cell{1};
a_raw = a_raw_cell{1};
cells_loc = cells_loc_cell{1}';
cells_drug = cells_drug_cell{1}';
cells_conc = cells_conc_cell{1}';
cells_slice= repmat(frame(1,1),length(b_cell{1}),1);
cells_ID   = cells_ID_cell{1}';
cells_intmass = cells_intmass_cell{1}';
cells_meanmass = cells_meanmass_cell{1}';
cells_finalmass = cells_finalmass_cell{1}';

% Second add other cells to the vector if you have multiple slices
if ExpFrames - SingleWindow +1>=2
    for i=2:ExpFrames - SingleWindow +1
        b = [b; b_cell{i}];
        b_raw = [b_raw; b_raw_cell{i}];
        a_raw = [a_raw; a_raw_cell{i}];
        cells_loc =[cells_loc; cells_loc_cell{i}'];
        cells_drug =[cells_drug; cells_drug_cell{i}'];
        cells_conc = [cells_conc; cells_conc_cell{i}'];
        cells_slice = [cells_slice; repmat(frame(1,i),length(b_cell{i}),1);];
        cells_ID    = [cells_ID; cells_ID_cell{i}'];
        %there is a fencepost error here?
        cells_intmass = [cells_intmass, cells_intmass_cell{i}'];
        cells_meanmass = [cells_meanmass, cells_meanmass_cell{i}'];
        cells_finalmass = [cells_finalmass, cells_finalmass_cell{i}'];
    end
end


[bOv,cells_locOv,cells_drugOv,cells_concOv] = ...
        sgrmaker_mean(tracks,tracksColHeaders,frame(1,1),frame(2,end),minmass,minp);
%% To correct old data run it from here
noC       = 2; %no of controls
conc_base = unique(cells_conc(cells_conc>0&cells_drug~=2)); %Return concentrations


%% Calculate mean data for the whole experiment for the same minp

b_meanOv    = zeros(max(cells_drugOv)-noC,length(conc_base));
b_medianOv  = zeros(max(cells_drugOv)-noC,length(conc_base));
b_stdOv     = zeros(max(cells_drugOv)-noC,length(conc_base));
drugOv      = zeros(max(cells_drugOv)-noC,length(conc_base));
concOv      = zeros(max(cells_drugOv)-noC,length(conc_base));
cellCountOv = zeros(max(cells_drugOv)-noC,length(conc_base));


for jj = noC + 1:max(cells_drugOv)
    for ii = 1:length(conc_base)
        
        x = cells_drugOv == jj & cells_concOv == conc_base(ii);
        if max(x)==0
            b_meanOv(jj-noC,ii) = nan;
            b_medianOv(jj-noC,ii) = nan;
            b_stdOv(jj-noC,ii)    = nan;
            drugOv(jj-noC,ii)     = nan;
            concOv(jj-noC,ii)     = nan;
            cellCountOv(jj-noC,ii)= nan;
        else
            b_meanOv(jj-noC,ii) = mean(bOv(x));
            b_medianOv(jj-noC,ii) = median(bOv(x));
            b_stdOv(jj-noC,ii)    = std(bOv(x));
            drugOv(jj-noC,ii)     = mean(cells_drugOv(x));
            concOv(jj-noC,ii)     = mean(cells_concOv(x));
            cellCountOv(jj-noC,ii)= length(bOv(x));
        end
    end
end

conc_control       = length(unique(cells_conc(cells_drug<=noC)));
b_mean_controlOv   = zeros(noC,conc_control);
b_median_controlOv = zeros(noC,conc_control);
b_std_controlOv    = zeros(noC,conc_control);
CountControlOv     = zeros(noC,conc_control);

for jj = 1:noC
    cont_conc = unique(cells_conc(cells_drug==jj));
    for zz = 1:length(cont_conc)
        x              = cells_drugOv== jj&cells_concOv==cont_conc(zz);
        if max(x)==0
            b_mean_controlOv(jj,zz)  = nan;
            b_median_controlOv(jj,zz)= nan;
            b_std_controlOv(jj,zz)   = nan;
            CountControlOv(jj,zz)    = nan;
        else
            b_mean_controlOv(jj,zz)  = mean(bOv(x));
            b_median_controlOv(jj,zz)= median(bOv(x));
            b_std_controlOv(jj,zz)   = std(bOv(x));
            CountControlOv(jj,zz)    = length(bOv(x));
        end
    end
end

b_mean_controlOv(b_mean_controlOv==0)  = nan;
b_median_controlOv(b_median_controlOv==0)= nan;
b_std_controlOv(b_std_controlOv==0)   = nan;
CountControlOv(CountControlOv==0)    = nan;


GR = zeros(size(b));
MinContConc1 = min(cells_concOv(cells_drugOv==1));
MinContConc2 = min(cells_concOv(cells_drugOv==2));
if ethanolcontrol > noC && ethanolcontrol <= max(cells_drug)
    GR = 2.^(b./b_mean_controlOv(2,cont_conc==MinContConc1))-1; % used to be b_mean_controlOv(2)
    GR(cells_drug==ethanolcontrol) = 2.^(b(cells_drug==ethanolcontrol)./b_mean_controlOv(1,cont_conc==MinContConc1))-1;
    GR(cells_drug==1) = 2.^(b(cells_drug==1)./b_mean_controlOv(1,cont_conc==MinContConc1))-1;
else
    GR = 2.^(b./b_mean_controlOv(2,cont_conc==MinContConc2))-1;
    GR(cells_drug==1) = 2.^(b(cells_drug==1)./b_mean_controlOv(1,cont_conc==MinContConc1))-1;
end


%% bin data together

binSlices = ceil(max(cells_slice)/binningRate);
cells_binSlices = zeros(size(cells_slice));
Oldtimes = mean(frame,1)/3;
times = zeros(binSlices,1);
for i = 1:binSlices
    x   = cells_slice <= binningRate*i & cells_slice > binningRate*(i-1);
    cells_binSlices(x) = i;
    if i == binSlices
        times(i) =mean([frame(1,binningRate*(i-1)+1);frame(2,end)])/3;
    else
        times(i) = mean([frame(1,binningRate*(i-1)+1);frame(2,binningRate*i)])/3;
    end
end

%%
j            = max(cells_drug)-noC;
i            = length(conc_base);
step         = 1e-4;
cont1_loc    = cells_drug == 1 & cells_conc == min(cells_conc(cells_drug==1));
cont2_loc    = cells_drug == 2 & cells_conc == min(cells_conc(cells_drug==2));
ran          = -0.1:step:.1;


GR_mean_control    = zeros(noC,length(cont_conc),binSlices);
b_mean_control     = zeros(noC,length(cont_conc),binSlices);
GR_median_control  = zeros(noC,length(cont_conc),binSlices);
b_median_control   = zeros(noC,length(cont_conc),binSlices);
GR_std_control     = zeros(noC,length(cont_conc),binSlices);
b_std_control      = zeros(noC,length(cont_conc),binSlices);
CountControl       = zeros(noC,length(cont_conc),binSlices);

for ss = 1:binSlices
    for jj = 1:noC
        for ii = 1:length(length(cont_conc))
        x              = cells_drug== jj & cells_binSlices ==ss& cells_conc==cont_conc(ii);
        GR_mean_control(jj,ii,ss) = mean(GR(x));
        b_mean_control(jj,ii,ss) = mean(b(x));
        GR_median_control(jj,ii,ss)= median(GR(x));
        b_median_control(jj,ii,ss)= median(b(x));
        GR_std_control(jj,ii,ss)   = std(GR(x));
        b_std_control(jj,ii,ss)   = std(b(x));
        CountControl(jj,ii,ss)    = length(b(x));
        ycont2          = pdf(fitdist(b(cont2_loc & cells_binSlices ==ss),'Kernel'),ran).*step;
        ycont1          = pdf(fitdist(b(cont1_loc & cells_binSlices ==ss),'Kernel'),ran).*step;
        hellcont(ss)    = 1/sqrt(2)*sqrt(sum((sqrt(ycont2)-sqrt(ycont1)).^2));
        end
    end
end

ToRThreshold = ceil(mean(hellcont)*10)/10;




%%
edges        = -0.1:0.001:0.15;
edgesGR      = 0:0.05:1;
j            = max(cells_drug)-noC;
i            = length(conc_base);
b_mean       = zeros(j,i,binSlices);
GR_mean      = zeros(j,i,binSlices);
b_median     = zeros(j,i,binSlices);
GR_median    = zeros(j,i,binSlices);
b_std        = zeros(j,i,binSlices);
GR_std       = zeros(j,i,binSlices);
drug         = zeros(j,i,binSlices);
conc         = zeros(j,i,binSlices);
cellCount    = zeros(j,i,binSlices);
hellinger    = zeros(j,i,binSlices);
hellingerGR  = zeros(j,i,binSlices);

for ss = 1:binSlices
    for jj = noC + 1:max(cells_drug)
    %     figure( 'Name', druglabels{jj} )
        for ii = 1:length(conc_base)
    %         subplot(3,2,ii)
            x = cells_drug == jj & cells_conc == conc_base(ii);
            y = x & cells_binSlices ==ss;
            if max(y)== 0
                b_mean(jj-noC,ii,ss) = nan;
                GR_mean(jj-noC,ii,ss) = nan;
                b_median(jj-noC,ii,ss) = nan;
                GR_median(jj-noC,ii,ss) = nan;

                b_std(jj-noC,ii,ss)    = nan;
                GR_std(jj-noC,ii,ss)    = nan;
                drug(jj-noC,ii,ss)     = nan;
                conc(jj-noC,ii,ss)     = nan;
                cellCount(jj-noC,ii,ss)= nan;
                hellinger(jj-noC,ii,ss)= nan;
                hellingerGR(jj-noC,ii,ss)= nan;
                
            else
                b_mean(jj-noC,ii,ss) = mean(b(y));
                GR_mean(jj-noC,ii,ss) = mean(GR(y));
                b_median(jj-noC,ii,ss) = median(b(y));
                GR_median(jj-noC,ii,ss) = median(GR(y));

                b_std(jj-noC,ii,ss)    = std(b(y));
                GR_std(jj-noC,ii,ss)    = std(GR(y));
                drug(jj-noC,ii,ss)     = mean(cells_drug(y));
                conc(jj-noC,ii,ss)     = mean(cells_conc(y));
                cellCount(jj-noC,ii,ss)= length(b(y));
                ydcb            = pdf(fitdist(b(y),'Kernel'),ran).*step;

                if jj~=ethanolcontrol
                    ycont2          = pdf(fitdist(b(cont2_loc & cells_binSlices ==ss),'Kernel'),ran).*step;
                    hellinger(jj-noC,ii,ss)= 1/sqrt(2)*sqrt(sum((sqrt(ydcb)-sqrt(ycont2)).^2));
                    
                else
                    ycont1          = pdf(fitdist(b(cont1_loc & cells_binSlices ==ss),'Kernel'),ran).*step;
                    hellinger(jj-noC,ii,ss)= 1/sqrt(2)*sqrt(sum((sqrt(ydcb)-sqrt(ycont1)).^2));

                end
            end
        end
    end
end


%%
polydegree   = 2;
hellpoly     = zeros(j,i,polydegree+1);
hellpolyplot = zeros(j,i,binSlices);
for jj = noC + 1:max(cells_drug)
    for ii = 1:length(conc_base)
        hellpoly(jj-noC,ii,:) =polyfit(times,squeeze(hellinger(jj-noC,ii,:)),polydegree);
        hellpolyplot(jj-noC,ii,:) = polyval(squeeze(hellpoly(jj-noC,ii,:)),times);
        if min(hellpolyplot(jj-noC,ii,:))>ToRThreshold
            ToR(jj-noC,ii)            = 0;
        elseif isnan(min(hellpolyplot(jj-noC,ii,:)))
            ToR(jj-noC,ii)            = nan;
        else
            ToR(jj-noC,ii)            = fsolve(@(t) squeeze(hellpoly(jj-noC,ii,3))-ToRThreshold+squeeze(hellpoly(jj-noC,ii,2))*t+squeeze(hellpoly(jj-noC,ii,1))*t^2,max(times));
        end
    end
end

%%
hellingerControl    = zeros(j,i,binSlices);

[sortedID, idx] = sort(cells_ID);
sortedSGR = b(idx);
sortedBins = cells_binSlices(idx);
sortedDrug = cells_drug(idx);
sortedConc = cells_conc(idx);
sortedLoc = cells_loc(idx);

clear idx

numBins = 18;

uniqueID = unique(sortedID);
meanSGR = zeros(length(uniqueID)*numBins,1);
abridgedBins = zeros(length(uniqueID)*numBins,1);
abridgedID = zeros(length(uniqueID)*numBins,1);
abridgedDrug = zeros(length(uniqueID)*numBins,1);
abridgedConc = zeros(length(uniqueID)*numBins,1);
abridgedLoc = zeros(length(uniqueID)*numBins,1);

counter = 1;
for ii = 1:length(uniqueID)
    currentDrug = sortedDrug(sortedID==uniqueID(ii));
    currentConc = sortedConc(sortedID==uniqueID(ii));
    currentLoc = sortedLoc(sortedID==uniqueID(ii));
    
    for jj = min(sortedBins(sortedID==uniqueID(ii))):max(sortedBins(sortedID==uniqueID(ii)))
        % finds mean SGR for a particular cellID and bin#
        meanSGR(counter) = mean(sortedSGR(sortedID==uniqueID(ii) & sortedBins==jj));
        abridgedBins(counter) = jj;
        abridgedID(counter) = uniqueID(ii);
        abridgedDrug(counter) = currentDrug(1);
        abridgedConc(counter) = currentConc(1);
        abridgedLoc(counter) = currentLoc(1);
        % below is unvalidated code, needs to be checked for accuracy
        counter = counter + 1;
    end
end

firstZero = find(abridgedLoc == 0, 1, 'first');

abridgedBins(firstZero:end) = [];
abridgedConc(firstZero:end) = [];
abridgedDrug(firstZero:end) = [];
abridgedID(firstZero:end) = [];
abridgedLoc(firstZero:end) = [];
meanSGR(firstZero:end) = [];
toc

% uncomment next line if you would like the results to be saved automatically
% save([froot, savename])
