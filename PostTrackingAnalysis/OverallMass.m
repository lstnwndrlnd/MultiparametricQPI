%%% This code 'OverallMass.m' is copyright 2021
% Authors: Tarek E. Mustafa, Edward R. Polanco, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

%%
close all, clear all
%%
% This code aims to  calculate Overall mass over time.

% froot should point to folder containing the raw data collected by the microscope
froot = '';
savename = froot;
panel = 2;
fstart= 'data';
fext  = '.mat';
filtersize = 5;     % Size of median filter element
ethanolcontrol = 0; % 3 if Panel A, 0 otherwise
wavelength = 624;   % Current wavelength
frames     = 187; % needs to be entered exactly!
% druglabels = {'EtOH','DMSO','4-HT','Lapt','Fulv','Doxo','5Fu'}; % For Panel A
druglabels = {'Untreated', 'DMSO', 'palbociclib','vinblastine','carboplatin','docetaxel',}; % for Panel B
%% Method 1 using D_stored and L_stored from data file (Most Raw form)

% 1) Load data file and calculate mass of the frame
wellnum       = 96;
Loc_mass      = zeros(frames,9*wellnum);
Loc_cellmass  =zeros(frames,9*wellnum);
Well_mass     = zeros(frames,wellnum);
Well_massstd  = zeros(frames,wellnum);


tic
for loc = 1:wellnum*9
    load([froot, fstart, num2str(loc),fext],'D_stored')

    if isa(D_stored,'single')
        Loc_mass(:,loc) = squeeze(sum(sum(double(D_stored(:,:,1:frames)).*wavelength,2),1));
    elseif isa(D_stored,'int16')
        Loc_mass(:,loc) = squeeze(sum(sum(double(D_stored(:,:,1:frames))*(2*pi).*wavelength/65536,2),1));
    else
        error('Error. \nInput must be a int16 or single, not a %s.',class(D_stored))
    end

    Loc_mass(:,loc) = medfilt1(Loc_mass(:,loc),filtersize,'truncate');     
end
    
    Loc_mass = Loc_mass(:,:)./Loc_mass(1,:);
    
for well =1:wellnum
    Well_mass(:,well) = mean(Loc_mass(:,9*well-8:9*well),2);
    Well_massstd(:,well)=std(Loc_mass(:,9*well-8:9*well),0,2);
end

load([froot,'data_allframes.mat'],'tracks','tracksColHeaders');
[~, cells_loc, cells_drug, cells_conc] = cellhunt_frame(tracks,[]);


cellsN_loc = histc(cells_loc,[1:864]);

cellsD_loc = zeros(1,864);
cellsC_loc = zeros(1,864);
cellsD_loc(1) = cells_drug(1);
cellsC_loc(1) = cells_conc(1);

for nn = 2:864
    if cellsN_loc(nn) ~= 0
        cellsD_loc(nn)=cells_drug(sum(cellsN_loc(1:nn-1))+1);
        cellsC_loc(nn)=cells_conc(sum(cellsN_loc(1:nn-1))+1);
        
    else
        if rem(nn-1,9) == 0
            cellsD_loc(nn)=cells_drug(sum(cellsN_loc(1:nn))+1);
            cellsC_loc(nn)=cells_conc(sum(cellsN_loc(1:nn))+1);
        else
            cellsD_loc(nn)=cellsD_loc(nn-1);
            cellsC_loc(nn)=cellsC_loc(nn-1);
        end
    end
end



for i = 1:96
    cellsD_well(i) = sum(cellsD_loc(9*i-8:9*i))/9;
    cellsC_well(i) = sum(cellsC_loc(9*i-8:9*i))/9;
    cellsN_well(i) = sum(cellsN_loc(9*i-8:9*i));
end


noC       = length(unique(cellsD_well(cellsC_well==0))); %no of controls
conc      = unique(cellsC_well(cellsC_well>0&cellsD_well~=2));

dcm = zeros(max(cellsD_well)-noC,length(conc),frames);
dcm_std = zeros(max(cellsD_well)-noC,length(conc),frames);
for jj = noC+1 : max(cellsD_well)
    for ii = 1 : length(conc)
        x  = cellsD_well == jj & cellsC_well == conc(ii);
        % dcm is drug % concentration mass
        if max(x) == 0
            
            dcm(jj-noC,ii,:) = nan;
            dcm_std(jj-noC,ii,:)   = nan;
        else
            dcm(jj-noC,ii,:) = mean(Well_mass(:,x),2);
            dcm_std(jj-noC,ii,:)   = std(Well_mass(:,x),0,2);
        end
    end
end

conc_control      = unique(cellsC_well(cellsD_well<=noC));

dcm_control = zeros(noC,length(conc_control),frames);
dcm_control_std = zeros(noC,length(conc_control),frames);

for jj = 1:noC
    for ii = 1:length(conc_control)
        x  = cellsD_well == jj & cellsC_well == conc_control(ii);
        if max(x) == 0
            dcm_control(jj,ii,:)  = mean(Well_mass(:,x),2);
            dcm_control_std(jj,ii,:) = std(Well_mass(:,x),0,2);
        else
            dcm_control(jj,ii,:)  = mean(Well_mass(:,x),2);
            dcm_control_std(jj,ii,:) = std(Well_mass(:,x),0,2);
        end
    end
end

times   = (1:frames)./3;

%%
save([savefolder, savename])
disp('workspace saved')
toc
datetime('now')
%%
for jj = noC+1 : max(cellsD_well)
    figure(jj)
    hold on
    if jj == ethanolcontrol
        errorbar(times,squeeze(dcm_control(1,1,:)),squeeze(dcm_control_std(1,1,:))./3)
    else
        errorbar(times,squeeze(dcm_control(2,1,:)),squeeze(dcm_control_std(2,1,:))./3)
    end
    for ii = 1 : length(conc)
        errorbar(times',squeeze(dcm(jj-noC,ii,:)),squeeze(dcm_std(jj-noC,ii,:))./3)
    end
    hold off
    ylim([0.8 7])
    title(druglabels{jj})
    if jj == ethanolcontrol
        labels = {druglabels{1}, num2str(conc(1)), num2str(conc(2)), ...
            num2str(conc(3)), num2str(conc(4)), num2str(conc(5)), ...
            num2str(conc(6))};
    else
        if panel ==1
            labels = {druglabels{2}, num2str(conc(1)), num2str(conc(2)), ...
                num2str(conc(3)), num2str(conc(4)), num2str(conc(5)), ...
                num2str(conc(6))};
        elseif panel == 2 && jj ~=6
                labels = {druglabels{2}, num2str(conc(end-5)), num2str(conc(end-4)), ...
                num2str(conc(end-3)), num2str(conc(end-2)), num2str(conc(end-1)), ...
                num2str(conc(end))};
        else
            labels = {druglabels{2}, num2str(conc(1)),num2str(conc(2)),...
                num2str(conc(end-4)), num2str(conc(end-3)), num2str(conc(end-2)),...
                num2str(conc(end-1)), num2str(conc(end))};
        end
    end
    legend(labels)
    xlabel('Time (h)', 'interpreter', 'latex' );
    ylabel( 'Normalized Mass (a.u.)', 'interpreter', 'latex');
end
