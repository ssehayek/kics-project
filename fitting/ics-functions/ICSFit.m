%
function out = ICSFit(params,xi_grid,eta_grid,varargin)

errBool = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        corr = varargin{i+1};
    end
end

% "getICSCorrSub.m" may replace (0,0) lag by NaN, so this value may need to
% be removed before performing least-squares
if errBool 
    nan_inds = find(isnan(corr));
    
    xi_grid(nan_inds) = [];
    eta_grid(nan_inds) = [];
    corr(nan_inds) = [];
end

s = struct('amp',params(1),'w0',params(2),'offset',params(3));
    
F = s.amp.*exp(-(xi_grid.^2+eta_grid.^2)/s.w0^2)+s.offset;

if errBool
    err = norm(corr - F,2);
end

if ~errBool
    out = F;
else
    out = err;
end
