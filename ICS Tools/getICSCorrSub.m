% This funciton returns a subset of the ICS function.
%
% INPUT PARAMS
% corr: entire ICS function with (0,0) lag centered using "fftshift.m"
% xi[eta]_lags: range of lags to include in the x[y]-direction
%               (specified as actual lag values e.g. 0 corresponds to
%               0th lag)
% 
% VARARGIN
% - choose whether to include (0,0) lag; value is set to NaN if true
%   'includeZeroLag' | values: (default 0) 1
function [sub_corr,xi_sub_grid,eta_sub_grid] = getICSCorrSub(corr,xi_lags,...
    eta_lags,varargin)

includeZeroLag = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'includeZeroLag','zeroLag'}))
        if varargin{i+1} == 0 || varargin{i+1} == 1
            includeZeroLag = varargin{i+1};
        else
            warning(['Unknown option for '' ',varargin{i},...
                ''', using default options.'])
        end
    end
end

% remove duplicates and sort lag vector inputs
xi_lags = unique(xi_lags);
eta_lags = unique(eta_lags);

% find (0,0) lag value position in "corr" input
[ctr_pxl_y,ctr_pxl_x] = getCtrPxl(corr);

[xi_sub_grid,eta_sub_grid] = meshgrid(xi_lags,eta_lags); % meshgrid subset
sub_corr = corr(ctr_pxl_y+eta_lags,ctr_pxl_x+xi_lags,:); % "corr" subset

if ~includeZeroLag % remove (0,0) lag and set it to NaN in lag grids and "sub_corr"
    [zero_i_y,zero_i_x] = find(xi_sub_grid == 0 & eta_sub_grid == 0);
    
    xi_sub_grid(zero_i_y,zero_i_x) = NaN;
    eta_sub_grid(zero_i_y,zero_i_x) = NaN;
    sub_corr(zero_i_y,zero_i_x,:) = NaN;
end