function [IC50, fitresult, gof] = hill4(concentrations, y_data, pcutoff)
% This function takes Concentration (x-data) and the response (y-data) to
% comupte IC50, fitresults, gof (check fit function help for more info) &
% figh which is figure handle. 
%This function uses 4 parameters Hill's equation to determine IC50.
% The function takes equal length array or vector of input
% To solve it for multiple conditions, put it in a loop

% Note that code will fail to compute EC50s because of d limits. TO compute
% solve either add a negative sign to d in the equation or flip to limits
% of the equation
% Input Limits and Starting points are found in sigfit equation and
% currently set to 
% a (max response depth, floor) Starting point is min of y_data and range
% is 5* Range of y_data
% b (min Response, ceiling) Starting point is max of y_data and the range
% is 5* Range of y_data
% c (IC50 concentration in -log10) Starting point is geomean of
% concentrations and Range is Range of concentrations
% d (hill slope) starting point is 1 and Range is from 0.1 to 10. Note: d is
% a downward slope because concentration is calculated in -log10. Min is
% set to 0.1 in order not to get a flat line
% Function uses an F-test to compare results to a flat fit to remove false
% fits. In case of failing to fit, IC50 is set to inf or -inf according to
% growth rate. 

%%% This code 'hill4.m' is copyright 2021
% Authors: Marc Hapfner, Mario Niepel, Mirra Chung, and Peter Sorger
% This code may be distributed in its original form when properly attributed.

% This code was originally intended to accompany the following articles 
% published in Nature Methods: 
% Alternative drug sensitivity metrics improve preclinical cancer pharmacogenomics
% doi: https://doi.org/10.1038/nbt.3882
% Growth rate inhibition metrics correct for confounders in measuring sensitivity to cancer drugs
% doi: https://doi.org/10.1038/nmeth.3853

% Modified November 2020 by Tarek E. Mustafa, Edward R. Polanco, and Thomas A.
% Zangle


if nargin ==2
    pcutoff = 0.05;
end

%% Fit: 'y_data'.
%prepareCurveData formats the data in a way that doesn't cause fit equation
%to fail. It's like a secertary preparing data for the boss
[xData, yData] = prepareCurveData( concentrations, y_data );

[fitresult_hill4, gof_hill4] = sig_fit( xData, yData );
[fitresult_flat , gof_flat ] = flat_fit( xData,yData );

% F-test for the significance of the sigmoidal fit
Npara_hill4 = 4; % N of parameters in the growth curve (4 parameters hill equation)
Npara_flat = 1; % Constant value of the flat line
RSS2 = gof_hill4.sse;
RSS1 = gof_flat.sse;
df1 = (Npara_hill4 -Npara_flat);
df2 = (length(yData) -Npara_hill4 +1);
F = ( (RSS1-RSS2)/df1 )/( RSS2/df2 );
pval = 1-fcdf(F, df1, df2);

if pval>=pcutoff || isnan(RSS2) %Failing the Significance test (flat fit)
    if fitresult_flat.a > 0.5
        IC50 = Inf;
    else
        IC50 = -Inf;
    end
    fitresult = fitresult_flat;
    
    gof       = gof_flat;
else %Significant fit (Sigmoidal fit)
    fitresult = fitresult_hill4;
    Parameters = coeffvalues(fitresult_hill4); 
    IC50 = 10^-Parameters(3);
    gof      = gof_hill4;
end

end

function [fitresult, gof] = sig_fit( concentrations, yData)
% Set up fittype and options.
% Parameters meaning: 
% a is the response at Dose = 0 & Starting point is the mean of y_data
% Limits for a are set for 5-folds the range of y_data
% b is the response at Dose = Inf & Starting point is the mean of y_data
% Limits for b are set for 5-folds the range of y_data
% c is the EC50 or IC50 according to the data & starting point is the
% geometric mean of concentrations
% Limits are one order of magnituide higher of max and one lower of min 
% concentration 
% d is the slope for the line transitioning from a to b.
aSP = min(yData);
aLimit = [min(yData)-range(yData)*2, max(yData)+range(yData)*2];
bSP = max(yData);
bLimit = [min(yData)-range(yData)*2, max(yData)+range(yData)*2];
cSP = -log10(median(concentrations));
cLimit = [min(concentrations)*10^-2 max(concentrations)*10^2];
cLimit = -log10(cLimit([2 1]));
dLimit = [0.1 5];
dSP = 1;
% 'a+(b-a)./
ft = fittype( 'a+(b-a)./(1+(x.*(10.^c)).^d)', 'independent', 'x', 'dependent', 'y' );

% Set up the method 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';

%limits and startpoint
opts.Lower = [aLimit(1) bLimit(1) cLimit(1) dLimit(1)];
opts.StartPoint = [aSP bSP cSP dSP];
opts.Upper = [aLimit(2) bLimit(2) cLimit(2) dLimit(2)];

% Fit model to data.
[fitresult, gof] = fit( concentrations, yData, ft, opts );

end

function [fit_result, gof] = flat_fit(concentrations, yData)
warning('off')
% flat fit function (special case of sigmoidal fit)
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
aSP = mean(yData);
aLimit = [min(yData)-range(yData)*2, max(yData)+range(yData)*2];
opts.Display = 'Off';
opts.Lower = [aLimit(1)];
opts.StartPoint = [aSP];
opts.Upper = [aLimit(2)];
f = fittype('a+0.*x', 'independent', 'x', 'dependent', 'y');
% Set up the method 
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
[fit_result, gof] = fit(concentrations, yData,f);
end

