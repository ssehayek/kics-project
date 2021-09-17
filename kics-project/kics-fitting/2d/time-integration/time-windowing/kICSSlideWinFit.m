% kICSSlideWinFit(...) help header

function out = kICSSlideWinFit(params,kSq,tauVector,T_s,varargin)

all_vars = {'D','r','K','frac','Ts','ksq','tau'};
all_vals = {params(1),params(2),params(3),params(4),T_s,kSq,...
    tauVector};

errBool = 0;
sym_vars = {''}; % initialize no variables as symbolic
for i = 1:length(varargin)
    % return output error (input calculated ACF)
    if any(strcmpi(varargin{i},{'error','err','residual','res'}))
        errBool = 1;
        ydata = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'symVars','symsVars'}))
        % choose symbolic variables
        if iscell(varargin{i+1}) && all(ischar([varargin{i+1}{:}]))
            sym_vars = varargin{i+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    end
end
s = setSymVars(all_vars,all_vals,sym_vars);

% s = struct('diffusion',params(1),'r',params(2),'K',params(3),...
%     'frac',params(4));

[tauGrid,kSqGrid] = meshgrid(tauVector,kSq);

% fit function computation
% [static_term,static_term_norm] = kICSSlideWinImmBleachFullFn(s,tauGrid);
% [diff_term,diff_term_norm] = kICSSlideWinDiffBleachFullFn(s,kSqGrid,tauGrid);
[static_term,diff_term,static_term_norm,diff_term_norm] = kICSSlideWinFn(s,kSqGrid,tauGrid);

% PSF in k-space
% Ik = exp(-s.w0^2.*kSqGrid/8);
% Ik0 = exp(-s.w0^2.*kSq(1)/8);

% normalized correlation function
% F = (Ik.^2.*(s.frac.*diff_term+(1-s.frac).*static_term))./(Ik0.^2.*(s.frac.*diff_term_norm+(1-s.frac).*static_term_norm));
F = (s.frac.*diff_term+(1-s.frac).*static_term)./(s.frac.*diff_term_norm+(1-s.frac).*static_term_norm);
% F = static_term./static_term_norm;
% F = diff_term./diff_term_norm;

% the best measure for the error, so far, seems to be to calculate the LS
% of each curve in tau individually, and then sum it. This is instead of
% calculating the norm point-by-point.
if errBool
%     err = sum(sqrt(sum((F-ydata).^2,1)),2);
    err = sum(sum((F-ydata).^2,1),2);
end

if ~errBool % output function
    out = double(F);
else % output error
    out = double(err);
end

if isnan(out)
    error(['Function is undefined for params: ',num2str(params),'.'])
end