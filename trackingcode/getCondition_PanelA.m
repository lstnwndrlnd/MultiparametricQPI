function [drugID, concentration, drugLabels, drugMap, concMap] = getCondition_PanelA(Loc)

% This code 'getCondition_PanelA.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

%define conditions in each well
drugidnums = [   1,      2,      3,      4,      5,       6,     7]';
% originally druglabels was a cell array (so you can have different length
% strings in each element). I changed it to an array, so only use four
% characters for each of the strings or this next line might break.
% drugLabels = ['EtOH'; 'DMSO'; '4-HT'; 'Lapt'; 'Fulv'; 'Doxo'; '5Fur'];
% drugLabels = {'ethanol', 'DMSO', '4-hydroxy-tamoxifen', 'fulvestrant', 'doxorubicin', '5-FU'};
drugLabels = {'Ethanol', 'DMSO', '4-Hydroxy-Tamoxifen','Lapatinib','Fulvestrant','Doxorubicin','5-Fluorouracil'};


% drugMap starts counting from 1. Once you find the row/col of the desired
% location, the corresponding row/col position has the drugID. The drugID
% is the index in drugLabels that translates drugID into a string.
drugMap_small = [1, 1, 1, 2, 2, 2, 3, 3, 3, 3, 3, 3;...
                 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3;...
                 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5;...
                 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5;...
                 4, 4, 4, 4, 4, 4, 5, 5, 5, 5, 5, 5;...
                 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7;...
                 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7;...
                 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7];

drugMap = imresize(drugMap_small, 3, 'nearest');

 conc_base = [.0016, .016, .08, .4, 2, 20]; % micromolar
 % repmat creates a 8x12 array of concentrations in the same pattern used
 % in the experiment.
 conc = repmat(conc_base,8,2);
 % in row 1, columns 1:6 are the controls. Set their concentration to 0.
 conc(1,1:6) = 0;
 
 concMap = imresize(conc, 3, 'nearest');

 % scale 8x12 array up to a 24x36 array. Use 'nearest method' to find what
 % value belongs in each element to just copy that element into neighboring
 % elements. This is because each of the 96 wells needs to be broken up
 % into 9 individual positions in the well. Since each position in each
 % well has the same condition, just copy the conditions into the empty
 % neighboring elements.
 % ex. [~ ~ ~    [ 1 1 1
 %      ~ 1 ~  =   1 1 1
 %      ~ ~ ~]     1 1 1 ]
 % the 3x3 empty matrix above is created by using imresize(~, 3, ~)
 % the 1 on the left represented an entire well before imresize()
 % the 3x3 of 1's on the right, represent each position in the first well
 % having the same condition as every other position in that well (in this
 % example the condition is drug 1, which is the etOH control)
 
%  drugs_all = imresize(drugid,3,'nearest');
%  conc_all = imresize(conc,3,'nearest');
 
 % Use help mapLocation or right click and open the function for a 
 % detailed explaination of how it works. Basically, it takes a location
 % numbered 1:864 and translates it into a corresponding row and column 
 % that can be used to find the specific condition used for that location.
 [locRow, locCol] = mapLocation(Loc);
 
 % once (row,col) coordinates are found, use them to identify condition.
 drugID = drugMap(locRow,locCol);
 concentration = concMap(locRow, locCol);


























end
