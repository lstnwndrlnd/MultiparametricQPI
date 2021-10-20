%%% This code 'IC50Calculator_DataAllFrames.m' is copyright 2021
% Authors: Tarek E. Mustafa, Edward R. Polanco, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

%% Clear workspace and closes windows
close all, clear all 

%% Inputs for the code

% Panel A:
druglabels = {'EtOH','DMSO','4-HT','Lapt','Fulv','Doxo','5Fu'}; % For Panel A,,kk,
panel = 1;

% New Panel B:
% druglabels = {'Untreated', 'DMSO', 'Palbociclib','Vinblastine','Carboplatin','Docetaxel'};
% panel = 2;

% froot should point to a folder containing data_allframes.mat
froot = './'; 

minp = 20;              % min path length for a track
% use StartTime and EndTime to decide which part of experiment to analyze (beginning, end, middle, etc)
% start and end time are in hours, and assumes 20 minutes in between frames
StartTime=0;            % time to start analyzing data (in hours) 
EndTime =72;            % take the portion of data you need (in hours)
growthThreshold = 0;    % specific Growth Rate Threshold
movaverage_points = 20; % calculate Overall Mass over a x points moving average
ethanolID         = 0;

if(panel == 1)
%     [drugID, conc, ~, ~, ~] = getCondition_PanelA_v1(Loc);
    druglabels = {'EtOH','DMSO','4-HT','Lapt','Fulv','Doxo','5Fu'}; % For Panel A
    ethanolID = 3;
elseif(panel == 2)
    druglabels = {'Untreated', 'DMSO', 'Palbociclib','Vinblastine','Carboplatin','Docetaxel'}; % for Panel B
%     [drugID, conc, ~, ~, ~] = getCondition_PanelB_v1(Loc);
end

% Drug ID whose control is ethanol
% use ethanolID to point to the condition using ethanol as a carrier 
% e.g. hydroxytamoxifen
% if there is no drug using ethanol as a carrier, point to one of the
% controls and the program will ignore this feature

% code will crash if we are not doing each experiment in
% triplicate on a plate

[cells, cells_loc, cells_drug, cells_conc,cells_ID] = cellhunt(froot, minp, StartTime, EndTime);
% This loop is used to preform regression for each cell that meet the min
% path length
leen = cellfun(@length,cells);

for i= 1:length(cells)
        % Set each cell inital time to zero
        cells_mod{i}= [cells{i}(:,1)-cells{i}(1,1) cells{i}(:,2)];
        
        % Calculate paramters used for regression
        % For more insight check Zar's book on biostatistics chapter 17
        xav(i)      = mean(cells_mod{i}(:,1));
        x            = cells_mod{i}(:,1)-xav(i);
        x2(i)       = sum(x.^2);
        yav(i)      = mean(cells_mod{i}(:,2));
        y            = cells_mod{i}(:,2)-yav(i);
        xy(i)       = sum(x.*y);
        b_raw(i)    = xy(i)/x2(i);
        % Growth rate of each cell
        a_raw(i)    = yav(i)-b_raw(i)*xav(i);
        % Y intercept of each cell (intial mass)
        r2(i)       = (xy(i)^2/x2(i))/sum(y.^2);
        TSS(i) = sum(y.^2);
        LRSS(i) = (xy(i)^2/x2(i));
        syx(i)  = sqrt((TSS(i) - LRSS(i))/(leen(i)-2));
        syxy(i) = syx(i)/yav(i);
        sb(i)   = sqrt(syx(i)^2/x2(i));
end

% normalize growth rate by y intercept
b           = b_raw./a_raw;
c=-1/(sqrt(2)*erfcinv(3/2));

%% Remove false cells

for jj = 1:max(cells_drug)
    drug_conc = unique(cells_conc(cells_drug==jj));
    for ii = 1:length(drug_conc)
        x = cells_drug == jj & cells_conc==drug_conc(ii) ;
        MADsyx = c*median(abs(syx(x)-median(syx(x))));
        MADb = c*median(abs(b(x)-median(b(x))));
        z = b>=median(b(x))+3*MADb | b<=median(b(x))-3*MADb;
        y = syx>=median(syx(x))+3*MADsyx;
        zy = z|y;
        b(x&zy) = [];
        cells_conc(x&zy)=[];
        cells_drug(x&zy)=[];
        cells_loc(x&zy)=[];
        cells_ID(x&zy)=[];
        syx(x&zy)=[];
        b_raw(x&zy)=[];
        a_raw(x&zy) =[];
        r2(x&zy)=[];
        cells(x&zy)=[];
        cells_mod(x&zy)=[];
    end
end
for i = 1:length(cells)
    rm(i) = cells{i}(end-2,2)/cells{i}(3,2);
    l2mass(i) = cells{i}(end-2,2);
    lastmass(i) = cells{i}(end,2);
    meanmass(i)     = mean(cells{i}(:,2));
end
leen = cellfun(@length,cells);

%% Reorganise data from tracks to be location specific and well specific
%Creates the cells arrays
%Legend for the cells arrays naming:
% N is number
% D is drug
% C is concentration
% Loc is per location 
% Well is per well
% e.g. 1- cellsN_loc is cells number per location (expect 864 entries)
% e.g. 2- cellsC_well is Concentration of drug per each well (expect 96
% entries)
cellsN_loc = histc(cells_loc,[1:864]);

b_vector =b;
cellsD_loc = zeros(1,864);
cellsC_loc = zeros(1,864);
cellsD_loc(1) = cells_drug(1);
cellsC_loc(1) = cells_conc(1);
% Marking each location Drug and concentration and checking if any
% location is empty
for nn = 2:864
    if cellsN_loc(nn) ~= 0
        cellsD_loc(nn)=unique(cells_drug((sum(cellsN_loc(1:nn-1))+1):sum(cellsN_loc(1:nn))));
        cellsC_loc(nn)=unique(cells_conc((sum(cellsN_loc(1:nn-1))+1):sum(cellsN_loc(1:nn))));
    else
        if rem(nn-1,9) == 0
            cellsD_loc(nn)=unique(cells_drug((sum(cellsN_loc(1:nn))+1):sum(cellsN_loc(1:nn+1)-1)));
            cellsC_loc(nn)=unique(cells_conc((sum(cellsN_loc(1:nn))+1):sum(cellsN_loc(1:nn+1)-1)));
        else
            cellsD_loc(nn)=cellsD_loc(nn-1);
            cellsC_loc(nn)=cellsC_loc(nn-1);
        end
    end
end


n = 1;
locs = 1:864;
for i = 1:96
    cellsD_well(i) = unique(cellsD_loc(9*i-8:9*i));
    cellsC_well(i) = unique(cellsC_loc(9*i-8:9*i));
    cellsN_well(i) = sum(cellsN_loc(9*i-8:9*i));
end

% Flipping code to match plate
 for ff = 1:6
   cellsD_well(ff*16-7:ff*16) = flip(cellsD_well(ff*16-7:ff*16));
   cellsC_well(ff*16-7:ff*16) = flip(cellsC_well(ff*16-7:ff*16));
   cellsN_well(ff*16-7:ff*16) = flip(cellsN_well(ff*16-7:ff*16));
 end
 
ncells = cellsN_well;
 % First iterations of the code had ncells as the array of number of cells
 % per well
 
 % Flipping code to match plate
 for ff = 1:6
   wellRange  = ff*16-7:ff*16;
   begin      = sum(ncells(1:ff*16-8))+1;
   ending     = sum(ncells(1:ff*16));
   b(begin:ending) = flip(b(begin:ending));
   cells_conc(begin:ending) = flip(cells_conc(begin:ending));
   cells_drug(begin:ending) = flip(cells_drug(begin:ending));
   cells_loc(begin:ending)  = flip(cells_loc(begin:ending));
%    ncells(ff*16-7:ff*16) = flip(ncells(ff*16-7:ff*16));
 end

for well= 1:96
    if well == 1
        begin      = 1;
        ending     = ncells(1);
        GrowthRate(well) = mean(b(begin:ending));
    else
        begin       = sum(ncells(1:well-1))+1;
        ending      = sum(ncells(1:well));
        GrowthRate(well) = mean(b(begin:ending));
    end
end

GrowthRate = reshape(GrowthRate,[8,12]);

drugid = reshape(cellsD_well,[8, 12]);
%Reshape the cells drug per well matrix to match each well with its drug
conc   = reshape(cellsC_well,[8, 12]);
cellNumber = reshape(cellsN_well,[8 12]);
conc_base = unique(conc(conc>0&drugid~=2));
drugidnums= unique(drugid);
noC       = 2; %no of controls


%Parameters for distribution graph
numIntervals    = 10;
maxSGR          = 0.05;
intervalWidth   = (2*maxSGR)/numIntervals;
x               = -0.04:intervalWidth:0.06;

% Parameters for Hellinger's distance
edges       = -0.1:0.001:0.15;

% SGR for controls
wells      = 1:96;
wells      = reshape(wells,[8,12]);

cc = [1 2];
for nn = 1:length(cc)
    cont_conc = unique(conc(drugid==nn));
    for zz = 1:length(cont_conc)
        loc          = wells(drugid==cc(nn)&conc==cont_conc(zz));
        b1           = b((sum(ncells(1:loc(1)-1))+1):(ncells(loc(1))+sum(ncells(1:loc(1)-1))));
        b2           = b((sum(ncells(1:loc(2)-1))+1):(ncells(loc(2))+sum(ncells(1:loc(2)-1))));
        b3           = b((sum(ncells(1:loc(3)-1))+1):(ncells(loc(3))+sum(ncells(1:loc(3)-1))));
        dcb          = [b1 b2 b3];
        b_control(cc(nn),zz*3-2:zz*3) = [mean(b1) mean(b2) mean(b3)];
        growing = dcb > growthThreshold;
        AboveLimit_control(cc(nn),zz) = length(find(growing))/length(dcb);
        ncount       = histc(dcb,x);
        RF_control(nn,zz,:)= ncount/length(dcb);
        std_control(nn,zz) = std(dcb);
        med_control(nn,zz) = median(dcb);
        sem_control(nn,zz) = std(dcb)/sqrt(length(dcb));
        ses_control(nn,zz)       = 1/sqrt(2*(length(dcb)-1));
        ncell_control(nn,zz)= length(dcb);
        dcb_control{nn,zz}  = dcb;
    end
end

%

hell(1) = 1/sqrt(2)*sqrt(sum((sqrt(histcounts(dcb,edges,'Normalization', 'probability'))-sqrt(histcounts(b1,edges,'Normalization', 'probability'))).^2));
hell(2) = 1/sqrt(2)*sqrt(sum((sqrt(histcounts(dcb,edges,'Normalization', 'probability'))-sqrt(histcounts(b2,edges,'Normalization', 'probability'))).^2));
hell(3) = 1/sqrt(2)*sqrt(sum((sqrt(histcounts(dcb,edges,'Normalization', 'probability'))-sqrt(histcounts(b3,edges,'Normalization', 'probability'))).^2));
hellb(1)=1/sqrt(2)*sqrt(sum((sqrt(histcounts(b1,edges,'Normalization', 'probability'))-sqrt(histcounts(b3,edges,'Normalization', 'probability'))).^2));
hellb(2)=1/sqrt(2)*sqrt(sum((sqrt(histcounts(b1,edges,'Normalization', 'probability'))-sqrt(histcounts(b2,edges,'Normalization', 'probability'))).^2));
hellb(3)=1/sqrt(2)*sqrt(sum((sqrt(histcounts(b2,edges,'Normalization', 'probability'))-sqrt(histcounts(b3,edges,'Normalization', 'probability'))).^2));
% plot(hell)
% hold on 
% plot(hellb)
% hold off

% SGR for each drug
bcont      = histcounts(dcb_control{2},edges,'Normalization', 'probability');
beth       = histcounts(dcb_control{1},edges,'Normalization', 'probability');

b_pop      = zeros(max(drugidnums)-max(noC),length(conc_base)*3);
med_pop    = zeros(max(drugidnums)-max(noC),length(conc_base));
std_pop    = zeros(max(drugidnums)-max(noC),length(conc_base));
sem_pop    = zeros(max(drugidnums)-max(noC),length(conc_base));
AboveLimit = zeros(max(drugidnums)-max(noC),length(conc_base));
hellinger  = zeros(max(drugidnums)-max(noC),length(conc_base));
% RF         = zeros(1,max(drugidnums)- max(noC),length(conc_base));
b_comb     = cell(max(drugidnums)-max(noC),length(conc_base));

% Rows for each drug
% column for each concentration
 for jj = 1:max(drugidnums)-noC

     for ii = 1:length(conc_base)
        loc   = [];
        loc   = wells(conc==conc_base(ii)&drugid==jj+noC);
        if isempty(loc)
                AboveLimit(jj,ii)  = nan;
                med_pop(jj,ii)     = nan;
                std_pop(jj,ii)     = nan;
                ncount             = nan;
%                 RF{jj}(ii,:)       = nan;
                ses_pop(jj,ii)     = nan;
                ncell_pop(jj,ii)   = nan;
                b_comb{jj,ii}      = nan;
                b_pop(jj,ii*3-2:ii*3) =  [nan nan nan];



                hellinger(jj,ii)   = nan;
        else
            if loc(1)==1
                b1   = b(1:ncells(1));
                b2   = b((sum(ncells(1:loc(2)-1))+1):(ncells(loc(2))+sum(ncells(1:loc(2)-1))));
                b3   = b((sum(ncells(1:loc(3)-1))+1):(ncells(loc(3))+sum(ncells(1:loc(3)-1))));
                dcb  = [b1 b2 b3];
            else
                b1   = b((sum(ncells(1:loc(1)-1))+1):(ncells(loc(1))+sum(ncells(1:loc(1)-1))));
                b2   = b((sum(ncells(1:loc(2)-1))+1):(ncells(loc(2))+sum(ncells(1:loc(2)-1))));
                b3   = b((sum(ncells(1:loc(3)-1))+1):(ncells(loc(3))+sum(ncells(1:loc(3)-1))));
                dcb  = [b1 b2 b3];
            end
                b_pop(jj,ii*3-2:ii*3) =  [mean(b1) mean(b2) mean(b3)];
                growing = dcb > growthThreshold;
                %Percentage of cells that above limit
                AboveLimit(jj,ii)  = length(find(growing))/length(dcb);
                med_pop(jj,ii)     = median(dcb);
                std_pop(jj,ii)     = std(dcb);
                ncount             = histc(dcb,x);
%                 RF{jj}(ii,:)       = ncount/length(dcb);
                ses_pop(jj,ii)     = 1/sqrt(2*(length(dcb)-1));
                ncell_pop(jj,ii)   = length(dcb);
                b_comb{jj,ii}       = dcb;

                hellinger(jj,ii)   = 1/sqrt(2)*sqrt(sum((sqrt(histcounts(dcb,edges,'Normalization', 'probability'))-sqrt(histcounts(b_comb{jj,1},edges,'Normalization', 'probability'))).^2));
            end
        end
 end

%% Calculate GR parameter for each drug and concentration parameter

if ethanolID > noC && ethanolID <= max(cells_drug)
    GRall = 2.^(b./trimmean(b(cells_drug==2&cells_conc==min(conc(drugid==2))),50))-1;
    GRall(cells_drug==ethanolID) = 2.^(b(cells_drug==ethanolID)./trimmean(b(cells_drug==1),50))-1;
    GRall(cells_drug==1) = 2.^(b(cells_drug==1)./trimmean(b(cells_drug==1),50))-1;
else
    GRall = 2.^(b./trimmean(b(cells_drug==2&cells_conc==min(conc(drugid==2))),50))-1;
    GRall(cells_drug==1) = 2.^(b(cells_drug==1)./trimmean(b(cells_drug==1),50))-1;
end
% GR = 2 ^ (SGR_Drug/SGR_control) - 1
xy      = size(b_pop);
GR      = zeros(xy(1),xy(2));

for jj = 1:max(drugidnums)-noC
    if jj == ethanolID
        controlno = 1;
    else
        controlno = 2;
    end
    for ii = 1:length(conc_base)
        GR(jj,ii*3-2:ii*3) = 2.^(b_pop(jj,ii*3-2:ii*3)/mean(b_control(controlno,b_control(controlno,:)~=0))) - 1; 
    end
end


%% IC50

EC50_SG       = zeros(max(drugidnums)-noC,1);
fitresults_SG = {};
EC50_GT       = zeros(max(drugidnums)-noC,1);
fitresults_GT = {};
GR50          = zeros(max(drugidnums)-noC,1);
fitresults_GR = {};
% EC50_RM       = zeros(max(drugidnums),1);
% fitresults_RM = {};

for i = 1:length(conc_base)
concrep(3*i-2:3*i) = conc_base(i);
end

for nn=1:max(drugidnums)-noC
    [EC50_SG(nn), fitresults_SG{nn}] = hill4(concrep(~isnan(b_pop(nn,:))),b_pop(nn,~isnan(b_pop(nn,:))),0.01);
    [EC50_GT(nn), fitresults_GT{nn}] = hill4(conc_base(~isnan(AboveLimit(nn,:)))',...
        AboveLimit(nn,~isnan(AboveLimit(nn,:))));
%     [EC50_RM(nn), fitresults_RM{nn}] = EC50(conc_base,RelMass(nn,:));
    [GR50(nn), fitresults_GR{nn}, DepthOfResponseGR(nn)]     = hill3(concrep(~isnan(GR(nn,:))),GR(nn,~isnan(GR(nn,:))),0.01);
end

DepthOfResponseGR = DepthOfResponseGR';
% rewrite b_pop
b_pop_rep = b_pop;

% uncomment next line to automatically save output
% save([froot, 'IC50_output_',num2str(StartTime),'_',num2str(EndTime)], '-v7.3')
toc
