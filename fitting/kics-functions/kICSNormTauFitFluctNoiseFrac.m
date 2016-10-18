% kICSNormTauFitFluctNoiseFrac(...) help header

function out = kICSNormTauFitFluctNoiseFrac(params,kSq,tauVector,varargin)

errBool = 0;
for i = 1:length(varargin)
    % return output error (input calculated ACF)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        ydata = varargin{i+1};
    end
end

s = struct('diffusion',params(1),'r',params(2),'K',params(3),...
    'frac',params(4),'w0',params(5),'sigma',params(6));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq);

A = s.diffusion.*kSqGrid;

% fit function computation
%
% diffusing autocorrelation factor
diff_term = (s.r+(1-s.r).*exp(-s.K.*tauGrid)).*exp(-A.*tauGrid);
%
% static autocorrelation term
static_term = (1-s.r).*exp(-s.K.*tauGrid);
%

% PSF in k-space
Ik = exp(-s.w0^2.*kSqGrid/8); 

% normalization
F_norm = Ik.^2.*(1-s.r.*(1-s.frac))+s.sigma;

% normalized correlation function
F = (Ik.^2.*(s.frac.*diff_term+(1-s.frac).*static_term))./...
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

if isnan(out)
    error(['Function is undefined for params: ',num2str(params),'.'])
end
