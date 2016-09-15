% kICSNormTauFitFluct(...) help header

function out = kICSNormTauFitFluct(params,kSq,tauVector,varargin)

errBool = 0;
tauNorm = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'})) % output error (input calculated ACF)
        errBool = 1;
        ydata = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'tauNorm','normByTau','normTau','normByLag','normLag'})) % choose tau normalization
        tauNorm = varargin{i+1};
    end
end

s = struct('diffusion',params(1),'k_off_frac',params(2),'K',params(3),'frac',params(4));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq); 

photophys_AC = 1-s.k_off_frac+s.k_off_frac*exp(-s.K*tauGrid); % photophysical factor multiplying diffusing correlation (state correlation)
photophys_norm = 1-s.k_off_frac+s.k_off_frac*exp(-s.K*tauNorm); % photophysical factor multiplying diffusing correlation in denom
        
% normalization
F_norm = s.frac.*photophys_norm.*exp(-kSqGrid.*s.diffusion*tauNorm)+(1-s.frac).*(photophys_norm-1);

% normalized correlation function
F = (s.frac*photophys_AC.*exp(-kSqGrid.*s.diffusion.*tauGrid)+(1-s.frac).*(photophys_AC-1))./...
    F_norm;

% the best measure for the error, so far, seems to be to calculate the LS
% of each curve in tau individually, and then sum it. This is instead of
% calculating the norm point-by-point.
if errBool
    err = sum(sqrt(sum((F-ydata(:,tauVector+1)).^2,1)),2);
end

if ~errBool % output function
    out = F;
else % output error
    out = err;
end