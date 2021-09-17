% timeIntkICSFit(...) help header

function out = timeIntkICSFit(params,kSq,tauVector,varargin)

errBool = 0;
for i = 1:length(varargin)
    % return output error (input calculated ACF)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        ydata = varargin{i+1};
    end
end

s = struct('diffusion',params(1),'r',params(2),'K',params(3),'frac',params(4));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq);

% fit function computation
[diff_term,static_term,diff_term_norm,static_term_norm] = ...
    timeIntFn(s,tauGrid,kSqGrid);

% normalization
F_norm = s.frac.*diff_term_norm+(1-s.frac).*static_term_norm;

% normalized correlation function
F = (s.frac.*diff_term+(1-s.frac).*static_term)./...
    F_norm;

% the best measure for the error, so far, seems to be to calculate the LS
% of each curve in tau individually, and then sum it. This is instead of
% calculating the norm point-by-point.
if errBool
    err = sqrt(sum(sum((F-ydata).^2,1),2));
end

if ~errBool % output function
    out = F;
else % output error
    out = err;
end

if isnan(out)
    error(['Function is undefined for params: ',num2str(params),'.'])
end