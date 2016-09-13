% kICSNormTauFitFluctNoiseBleach(...) help header

function out = kICSNormTauFitFluctNoiseBleach(params,kSq,tauVector,k_p,T,varargin)

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

% photophysical factor multiplying diffusing correlation
photophys_AC = K^2/s.k_on^2*exp(1).^(k_p+(-1).*k_p.*T+(-1).*k_p.*tauGrid).*((-1)+exp(1).^k_p).^(-1).*( ...
    exp(1).^(k_p.*T)+(-1).*exp(1).^(k_p.*tauGrid)).*s.k_on.*K.^(-2).*( ...
    exp(1).^(((-1).*s.k_off+(-1).*s.k_on).*tauGrid).*s.k_off+s.k_on).*(T+(-1).*tauGrid).^(-1);
%
% photophysical factor multiplying static correlation
photophys_fluct_AC = K^2/s.k_on^2*(1/2).*exp(1).^(k_p+(-2).*k_p.*T+(-1).*(K+k_p).*tauGrid).*(exp(1).^( ...
    k_p.*T)+(-1).*exp(1).^(k_p.*tauGrid)).*s.k_on.*K.^(-2).*((-1).*exp( ...
    1).^(k_p+(K+k_p).*tauGrid).*s.k_on+exp(1).^(k_p.*T).*(s.k_off+exp(1).^k_p.* ...
    s.k_off+exp(1).^(K.*tauGrid).*s.k_on)).*(T+(-1).*tauGrid).^(-1).*((-1)+ ...
    coth(k_p));
%
% photophysical factor multiplying diffusing correlation in norm factor
photophys_AC_norm = K^2/s.k_on^2*exp(1).^(k_p+(-1).*k_p.*T+(-1).*k_p.*tauNorm).*((-1)+exp(1).^k_p).^(-1).*( ...
    exp(1).^(k_p.*T)+(-1).*exp(1).^(k_p.*tauNorm)).*s.k_on.*K.^(-2).*( ...
    exp(1).^(((-1).*s.k_off+(-1).*s.k_on).*tauNorm).*s.k_off+s.k_on).*(T+(-1).*tauNorm).^(-1);
%
% photophysical factor multiplying static correlation in norm factor
photophys_fluct_AC_norm = K^2/s.k_on^2*(1/2).*exp(1).^(k_p+(-2).*k_p.*T+(-1).*K.*tauNorm).*(exp(1).^( ...
    k_p.*T)+(-1).*exp(1).^(k_p.*tauNorm)).*s.k_on.*K.^(-2).*((-1).*exp( ...
    1).^(k_p+(K+k_p).*tauNorm).*s.k_on+exp(1).^(k_p.*T).*(s.k_off+exp(1).^k_p.* ...
    s.k_off+exp(1).^(K.*tauNorm).*s.k_on)).*(T+(-1).*tauNorm).^(-1).*((-1)+ ...
    coth(k_p));
%

Ik = exp(-s.w0^2/8*kSqGrid); % PSF in k-space

% normalization
F_norm = Ik.^2.*(s.frac.*photophys_AC_norm.*exp(-kSqGrid.*s.diffusion*tauNorm)+(1-s.frac).*photophys_fluct_AC_norm)+s.sigma;

% normalized correlation function
F = Ik.^2.*(s.frac*photophys_AC.*exp(-kSqGrid.*s.diffusion.*tauGrid)+(1-s.frac).*photophys_fluct_AC)./...
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