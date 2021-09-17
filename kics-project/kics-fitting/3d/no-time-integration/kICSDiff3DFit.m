% kICSDiff3DFit(...) help header

function out = kICSDiff3DFit(params,kSq,tauVector,varargin)

errBool = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'})) % output error (input calculated ACF)
        errBool = 1;
        ydata = varargin{i+1};
    end
end

s = struct('D',params(1),'r',params(2),'K',params(3),'frac',params(4),'z0',params(5));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq);

cxn_3d = 1./sqrt(4.*s.D.*tauGrid+s.z0.^2);

% normalization
F_norm = (s.frac+(1-s.frac).*(1-s.r))./s.z0;

% normalized correlation function
F = (s.frac.*(s.r+(1-s.r).*exp(-s.K.*tauGrid)).*exp(-kSqGrid.*s.D.*tauGrid).*cxn_3d+(1-s.frac).*(1-s.r).*exp(-s.K.*tauGrid)./s.z0)./...
    F_norm;


% the best measure for the error, so far, seems to be to calculate the LS
% of each curve in tau individually, and then sum it. This is instead of
% calculating the norm point-by-point.
if errBool
    err = sum(sqrt(sum((F-ydata).^2,1)),2);
end

if ~errBool % output function
    out = F;
else % output error
    out = err;
end