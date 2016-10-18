% timeIntkICSNoBleachFit(...) help header

function out = timeIntkICSFitNoBleach(params,kSq,tauVector,T,varargin)

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

F = exp(1).^((-1).*s.K.*(1+tauGrid)+(-1/4).*kSqGrid.*s.w0.^2).*(2.*exp(1).^((-1).*s.K+( ...
  -1/4).*kSqGrid.*s.w0.^2).*(((-1)+s.frac).*(1+exp(1).^s.K.*((-1)+s.K)).*s.K.^(-2).*(( ...
  -1)+s.r)+A.^(-2).*exp(1).^((-1).*A).*s.frac.*(A+s.K).^(-2).*(A.^3.*exp(1) ...
  .^(A+s.K)+A.*exp(1).^s.K.*(2+exp(1).^A.*((-2)+s.K)).*s.K.*s.r+(-1).*exp(1) ...
  .^s.K.*((-1)+exp(1).^A).*s.K.^2.*s.r+A.^2.*(1+(-1).*s.r+exp(1).^s.K.*s.r+exp( ...
  1).^(A+s.K).*((-1)+s.K+s.K.*s.r))))+s.sigma).^(-1).*(((-1)+exp(1).^s.K).^2.*((-1)+ ...
  s.frac).*s.K.^(-2).*((-1)+s.r)+(-1).*A.^(-2).*exp(1).^((-1).*A.*(1+tauGrid)).*s.frac.* ...
  (A+s.K).^(-2).*((-1).*A.^2.*((-1)+exp(1).^(A+s.K)).^2.*((-1)+s.r)+exp(1) ...
  .^(s.K+s.K.*tauGrid).*((-1)+exp(1).^A).^2.*(A+s.K).^2.*s.r).*(T+(-1).*tauGrid).^(-1).* ...
  ((-1).*T+tauGrid));

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
    error(['Function is undefined for params: ',num2str(params)])
end
