function [hellcont, Newhellinger] = HellingerDistComputation(abridgedDrug, abridgedConc, abridgedBins, ...
                                    meanSGR, conc_base, noC, binSlices, ethanolcontrol)
%%% This code 'HellingerDistComputation.m' is copyright 2021
% Authors: Tarek E. Mustafa, Edward R. Polanco, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

% tic
j            = max(abridgedDrug)-noC;
i            = length(conc_base);
step         = 1e-4;
cont1_loc    = abridgedDrug == 1 & abridgedConc == min(abridgedConc(abridgedDrug==1));
cont2_loc    = abridgedDrug == 2 & abridgedConc == min(abridgedConc(abridgedDrug==2));
ran          = -0.1:step:.1;

Newhellinger    = zeros(j,i,binSlices);
hellcont        = zeros(1,binSlices);
for ss = 1:binSlices
    for jj = noC + 1:max(abridgedDrug)
        for ii = 1:length(conc_base)
            x = abridgedDrug == jj & abridgedConc == conc_base(ii);
            y = x & abridgedBins ==ss;
            if max(y)==0 %isempty(y) | 

                Newhellinger(jj-noC,ii,ss)= nan;
            else
                ydcb            = pdf(fitdist(meanSGR(y),'Kernel'),ran).*step;
                ycont2          = pdf(fitdist(meanSGR(cont2_loc & abridgedBins ==ss),'Kernel'),ran).*step;
                ycont1          = pdf(fitdist(meanSGR(cont1_loc & abridgedBins ==ss),'Kernel'),ran).*step;
                hellcont(ss)    = 1/sqrt(2)*sqrt(sum((sqrt(ycont2)-sqrt(ycont1)).^2));
                if jj==ethanolcontrol
                    Newhellinger(jj-noC,ii,ss)= 1/sqrt(2)*sqrt(sum((sqrt(ydcb)-sqrt(ycont1)).^2));
                else
                    Newhellinger(jj-noC,ii,ss)= 1/sqrt(2)*sqrt(sum((sqrt(ydcb)-sqrt(ycont2)).^2));
                end
            end
        end
    end
end


% toc