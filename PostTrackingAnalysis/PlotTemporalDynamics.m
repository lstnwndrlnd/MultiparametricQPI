% This code 'PlotTemporalDynamics.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

% froot should point to output from sliding window code
froot = '.\slidingwindow_MDA_PanelB.mat';
load(froot);
%%
drugOfInterest = 6;
concOfInterest = 4;

meanControl = zeros(max(abridgedBins),1);
meanDrug = zeros(max(abridgedBins),5);

for ii = 1:max(abridgedBins)    
    meanControl(ii,1)    = nanmedian(meanSGR(abridgedDrug == 1 & abridgedConc == 0 & abridgedBins == ii));
end

for ii = 1:max(abridgedBins)
    meanDrug(ii,1)    = nanmedian(meanSGR(abridgedDrug == drugOfInterest & abridgedConc == conc_base(concOfInterest) & abridgedBins == ii));
end


figure
vl = violinplotleft(meanSGR(abridgedDrug == 2 & abridgedConc == 0), ...
    abridgedBins(abridgedDrug == 2 & abridgedConc == 0),            ...
    'ViolinColor',[0.3010 0.7450 0.9330]);

vr = ...
    violinplotright(meanSGR(abridgedDrug == drugOfInterest & abridgedConc == conc_base(concOfInterest)), ...
    abridgedBins(abridgedDrug == drugOfInterest & abridgedConc == conc_base(concOfInterest)),   ...
    'ViolinColor',[162, 29,  49]./255);    

lc = plot([1:max(abridgedBins)]'-0.05,meanControl(:,1),'Color',[0.3010 0.7450 0.9330]);
lc.LineWidth = 2;
lht = plot([1:max(abridgedBins)]'+0.03,meanDrug(:,1),'Color',[162, 29,  49]./255);
lht.LineWidth = 2;
xticklabels(round(times,1))
xlabel('Time (h)')
ylabel( 'Specific Growth Rate $(1/h)$', 'interpreter', 'latex');
xlabel('Time $(h)$','interpreter','latex')
title(['MDA-MB-231, control (left) and ', druglabels{drugOfInterest}, ' ', num2str(conc_base(concOfInterest)), '\muM (right) population over time'])

ylim([-.07 .12])