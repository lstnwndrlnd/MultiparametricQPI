% function [num, indices, lengths] = track_numpart(tracks, minL)
function [num, indices] = track_numpart(tracks, minL)
% returns the number of particles in the coordinate vector tracks. optional:
% minL specifies the min path length to count (ex. if minL = 2, then this
% function will only return the number of particles with x, y, z, t info at
% at least two time points
% inputs: tracks, array containing x, y, z, t, i (locations x, y, third
% coordinate z, and time coordinate, t for i particles); minL, optional
% input argument specifying the minimum path length to count
% outputs: num, the number of particles counted; indices, array of first 
% indices of %particles with path length > minL (if no minL specified, then
% it is automatically set to 0 to count all paths); lengths - path lengths
% (number of time points) for each particle labeled in indices

% This code 'track_numpart.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

if nargin == 1
    minL = 0;
end


lengths = [];
indices = [];
num = 0;
ind_i = 1;
ii = 1;
while ii <= length(tracks(:,1))
    len = track_lengthn(tracks, ii);
    if len >= minL
        num = num+1;
        lengths(ind_i) = len;
        indices(ind_i) = ii;
        ind_i = ind_i+1;
        ii = ii+len;
    else
        ii = ii+1;
    end
end


function [len] = track_lengthn(tracks, n)
%finds the length of the particle starting at location n
len = length(find(tracks(:,5)==tracks(n,5)));