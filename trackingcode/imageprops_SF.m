function [V, M, A, MI, P, SF] = imageprops_SF(L, I, pxlsize)
% function [V, M, A, MI, P, SF] = imageprops_SF(L, I, pxlsize)
% function to return volume, mass and area of regions L in image I
% inputs: L, the label image *can also be BW mask image. should be nonzero
% in regions where the image will be analyzed; I, the image to be processed;
% pxlsize, the size of each pixel in the image (mm/pixel)
% outputs: V, the measured volume (um^3) of each region; M the measured mass
% (pg); A, the area in pixels of each region; MI, the mean intensity of each
% region; P, struct containing all data returned by the call to the
% regionprops function ('Area', 'Centroid', 'MeanIntensity', 'BoundingBox',
% 'PixelValues') and the standard deviation ('Std')
% _SF version also returns 'shape factor' defined as the circularity
% (= 4pi*A/P^2)

% This code 'imageprops_SF.m' is copyright 2021
% Authors: Edward R. Polanco, Tarek E. Mustafa, Thomas A. Zangle
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

% TAZ modified 8/26/20 to remove imkeeplargest (this is slow)
% this may cause some issues with shape factor calculation
% also adjusted conditions for computing shape factor to avoid inf values
% and avoid calc for small regions (smaller than 9 pixels)

% ERP modified 1/20/21 to adapt code to convert from radians to phase shift
% in um for DPC microscope

%%% define relationship between mass and "volume"
K = 1./(10000).^3./100./0.0018.*1e12/(pi); %pg/um^3

P = regionprops(L, I, 'Area', 'Centroid', 'MeanIntensity', 'BoundingBox', 'PixelValues', 'Perimeter');

if max(max(L))>0
    for ii = 1:length(P)
        A(ii) = P(ii).Area;
        MI(ii) = P(ii).MeanIntensity;
        P(ii).Std = std(P(ii).PixelValues);
        if ~isempty(P(ii)) && P(ii).Perimeter>0 && P(ii).Area>9 %only compute if label exists, perimeter is nonzero and area is great enough for the calc to make sense
            SF(ii) = P(ii).Area*4*pi/(P(ii).Perimeter.^2);
        else
            SF(ii) = NaN;
        end
    end
else
    A = [];
    MI = [];
    P = [];
    SF = [];
end

V = MI.*A.*pxlsize.^2.*1e3/pi;

M = V.*K;
