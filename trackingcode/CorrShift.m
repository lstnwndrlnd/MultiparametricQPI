function [ishift, jshift] = CorrShift(A,B)
% find correlation shift between two matrices of the same size, A and B
% [ishift,jshift] is the amount to move elements of A such that they line up
% with elements of B

% This code 'CorrrShift.m' is copyright 2021
% Authors: Thomas A. Zangle, Edward R. Polanco, Tarek E. Mustafa
% This code is distributed under a creative commons attributable
% sharealike license. This license allows you to remix, adapt, and build 
% upon this work, as long as the authors are credited and the modified code
% is redistributed under the same license.

A(isnan(A)) = 0;
B(isnan(B)) = 0;

C = normxcorr2(A,B);

idim = length(A(:,1));
jdim = length(A(1,:));

[YY,II] = max(C);
[Y,J] = max(YY);

ishift = II(J)-idim;
jshift = J-jdim;