% find coordinates of (0,0) lag value in "corr" input, where
% "corr" is a correlation function which has its (0,0) lag at the centre
% after using "fftshift.m" on it. "corr" can include time (i.e. dimension
% of 3), since only the first two dimensions are considered in this
% function.
%
function [ctr_pxl_y,ctr_pxl_x] = getCtrPxl(corr)

size_y = size(corr,1);
size_x = size(corr,2);

% treat even and odd cases separately
if ~mod(size_y,2) % even
    ctr_pxl_y = size_y/2+1;
else % odd
    ctr_pxl_y = (size_y+1)/2;
end

if ~mod(size_x,2)
    ctr_pxl_x = size_x/2+1;
else
    ctr_pxl_x = (size_x+1)/2;
end
