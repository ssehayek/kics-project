% kICSNormTauFitFluctNoise(...) help header

function out = kICSNormTauFitFluctNoise(params,kSq,tauVector,varargin)

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

s = struct('diffusion',params(1),'k_on',params(2),'k_off',params(3),'frac',params(4),'w0',params(5),'sigma',params(6));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq); 

K = s.k_on + s.k_off; % sum of photophysical rates

photophys_AC = 1+s.k_off/s.k_on*exp(-K*tauGrid); % photophysical factor multiplying diffusing correlation (state correlation)
photophys_norm = 1+s.k_off/s.k_on*exp(-K*tauNorm); % photophysical factor multiplying diffusing correlation in denom

Ik = exp(-s.w0^2/8*kSqGrid); % PSF in k-space
        
% normalization
F_norm = Ik.^2.*(s.frac.*photophys_norm.*exp(-kSqGrid.*s.diffusion*tauNorm)+(1-s.frac).*(photophys_norm-1))+s.sigma;

% normalized correlation function
F = Ik.^2.*(s.frac*photophys_AC.*exp(-kSqGrid.*s.diffusion.*tauGrid)+(1-s.frac).*(photophys_AC-1))./...
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