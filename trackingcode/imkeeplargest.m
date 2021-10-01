function [ Lout ] = imkeeplargest( Lin )
%function [ Lout ] = imkeeplargest( Lin )
%function to find the only keep the largest regions in a label image, Lin
%Lin and Lout should be "label" images, with discrete, integer values
%starting at 0 (background) up to max(Lin(:)), the last region in the image

% This code 'imkeeplargest.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

Lout = zeros(size(Lin));

for ii = 1:max(Lin(:))
    Ltemp = bwlabel(Lin==ii);
    if sum(Ltemp(:)) %only continue if there is anything with this given label index
        P = regionprops(Ltemp, 'Area');
        [~,MaxIdx] = max(struct2array(P));
        Lout(Ltemp==MaxIdx) = ii;
    end


end

