%%% This code 'ComputeHellDistAndToRatEC50.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

experiment = 'C:\Users\Hassan\Documents\Data\DoDDataAnalysis\MDAMB231\PanelB\030221\slidingwindow_MDA_PanelB_73_8_030221.mat';
load(experiment)

%%
[hellcont Newhellinger] = HellCorrect_abridgedData(abridgedDrug, abridgedConc, abridgedBins, ...
                                    meanSGR, conc_base, noC, binSlices, ethanolcontrol);
%%

colorC = [255, 215,  0 ;    ... % yellow - 4HT
          127,  0,  255;    ... % purple - lap
           0 , 255,  0 ;    ... % green - fulv
          255,  0 ,  0 ;    ... % red - dox 
           89,  89,  91;    ... % blue - 5FU
           0 , 255, 255;    ... % cyan - palbo
          255,  0,  255;    ... % magenta - vinblastine
          255, 127,  0 ;    ... % orange - carboplatin
          162, 29,  49 ;    ... % scarlet - docetaxel
         ]./255;
     
figure
hold on
for ii = [2, 4]
    if ii == 2
        concOfInterest = 9; % EC50 concentration for vinblastine
    elseif ii == 4
        concOfInterest = 4; % EC50 concentration for docetaxel
    end
    
    scatter(times, squeeze(Newhellinger(ii,concOfInterest,:)), 100, ...
        'markerfacecolor', colorC(ii+5,:), 'markeredgecolor', 'none')
    s = squeeze(Newhellinger(ii,concOfInterest,:));

    fun = @(x,t)x(1)-x(2)*exp(x(3).*times); % model to fit to data
    x0 = [-.3379 -.5151 .01]; % initial guess at fitting weights
    lb = [0 0 -1]; % lower boundary constraint fitting weights
    ub = [1 1 0]; % upper boundary constraint fitting weights
    x = lsqcurvefit(fun, x0, times, s, lb, ub);
    
    xdata = linspace(0, 100, 500);
    plot(xdata, x(1)-x(2)*exp(x(3).*xdata), 'color', colorC(ii+5,:))

end
axis square
title('MDA-MB-231 treated with vinblastine (purple) and docetaxel (scarlet)')

hold off