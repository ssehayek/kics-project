function out = testTimeIntkICSDiffOnly(params,ksq,tau_vec,varargin)

errBool = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'})) % output error (input calculated ACF)
        errBool = 1;
        ydata = varargin{i+1};
    end
end
    
s = struct('D',params(1));

[tau_grid,ksq_grid] = meshgrid(tau_vec,ksq);

R_tau = exp(-s.D.*ksq_grid.*tau_grid).*(exp(s.D.*ksq_grid)-1).^2;
R_0 = 2+2*exp(s.D.*ksq_grid).*(s.D.*ksq_grid-1);

F = R_tau./R_0;

% the best measure for the error, so far, seems to be to calculate the LS
% of each curve in tau individually, and then sum it. This is instead of
% calculating the norm point-by-point.
if ~errBool % output function
    out = F;
else % output error
    err = sum(sqrt(sum((F-ydata).^2,1)),2);
    out = err;
end