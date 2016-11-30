%
function out = ICSFitAiry(params,xi_grid,eta_grid,varargin)

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
    
I_airy = @(x,y) (2.*besselj(1,s.w0.*sqrt(x.^2+y.^2))./sqrt(x.^2+y.^2)).^2;
% I_mult = @(ux,uy) I_airy(x-ux,y-uy).*I_airy(x+xi_grid-ux,y+eta_grid-uy);
I_mult = @(ux,uy,xi,eta) I_airy(ux,uy).*I_airy(xi-ux,eta-uy);

I_corr = zeros(size(xi_grid));
for ii = 1:numel(xi_grid)
    I_corr(ii) = integral2(@(ux,uy) I_mult(ux,uy,xi_grid(ii),eta_grid(ii)),...
        -100,100,-100,100);
end
I_ACF = I_corr/(4*pi);

F = s.amp.*I_ACF+s.offset;

if errBool
    err = norm(corr - F,2);
end

if ~errBool
    out = F;
else
    out = err;
end
