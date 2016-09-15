% kICSNormTauFitFluctNoiseBleach2(...) help header

function out = kICSNormTauFitFluctNoiseBleachBlinkFrac(params,kSq,tauVector,kp,T,varargin)

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

s = struct('diffusion',params(1),'r',params(2),'K',params(3),'frac',params(4),'w0',params(5),'sigma',params(6));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq);

% photophysical factor multiplying diffusing correlation
photophys_AC = exp(1).^((-1).*s.K.*tauGrid+(-1).*kp.*((-1)+T+tauGrid)).*((-1)+exp(1).^kp).^(-1) ...
  .*(exp(1).^(kp.*T)+(-1).*exp(1).^(kp.*tauGrid)).*(1+((-1)+exp(1).^( ...
  s.K.*tauGrid)).*s.r).*(T+(-1).*tauGrid).^(-1);
%
% photophysical factor multiplying static correlation
photophys_fluct_AC = (1/2).*exp(1).^((-2).*kp.*T+(-1).*(s.K+kp).*tauGrid).*(exp(1).^(kp.*T)+( ...
  -1).*exp(1).^(kp.*tauGrid)).*(exp(1).^(kp.*T).*(1+exp(1).^kp).*((-1)+ ...
  s.r)+(-1).*exp(1).^(s.K.*tauGrid).*(exp(1).^(kp.*T)+(-1).*exp(1).^(kp+kp.*tauGrid) ...
  ).*s.r).*((-1).*T+tauGrid).^(-1).*csch(kp);
%
% photophysical factor multiplying diffusing correlation in norm factor
photophys_AC_norm = exp(1).^((-1).*s.K.*tauNorm+(-1).*kp.*((-1)+T+tauNorm)).*((-1)+exp(1).^kp).^(-1) ...
  .*(exp(1).^(kp.*T)+(-1).*exp(1).^(kp.*tauNorm)).*(1+((-1)+exp(1).^( ...
  s.K.*tauNorm)).*s.r).*(T+(-1).*tauNorm).^(-1);
%
% photophysical factor multiplying static correlation in norm factor
photophys_fluct_AC_norm = (1/2).*exp(1).^((-2).*kp.*T+(-1).*(s.K+kp).*tauNorm).*(exp(1).^(kp.*T)+( ...
  -1).*exp(1).^(kp.*tauNorm)).*(exp(1).^(kp.*T).*(1+exp(1).^kp).*((-1)+ ...
  s.r)+(-1).*exp(1).^(s.K.*tauNorm).*(exp(1).^(kp.*T)+(-1).*exp(1).^(kp+kp.*tauNorm) ...
  ).*s.r).*((-1).*T+tauNorm).^(-1).*csch(kp);
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
    err = sum(sqrt(sum((F-ydata(:,tauVector+1)).^2,1)),2);
end

if ~errBool % output function
    out = F;
else % output error
    out = err;
end