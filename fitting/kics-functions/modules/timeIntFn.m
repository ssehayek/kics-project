function [diff_term,static_term,diff_term_norm,static_term_norm] = ...
    timeIntFn(s,tauGrid,kSqGrid)

A = s.diffusion.*kSqGrid;

% diffusing autocorrelation factor
diff_term = exp(-(s.K+A).*(1+tauGrid)).*(exp(s.K.*(1+tauGrid)).*(exp(A)-1).^2.*s.r./A.^2 + ...
    (1-s.r)./(A+s.K).^2.*(exp(A+s.K)-1).^2);
%
% static autocorrelation term
static_term = (1-s.r)./s.K.^2.*exp(-s.K.*(1+tauGrid)).*(exp(s.K)-1).^2;
%
% diffusing autocorrelation term (norm)
diff_term_norm = 2.*exp(-(s.K+A)).*1./(A.^2.*(A+s.K).^2).*(A.^2.*(1-s.r) + ...
    exp(s.K).*((A+s.K).^2.*s.r+exp(A).*(A.^2.*(A+s.K-1)+s.K.*s.r.* ...
    (-s.K+A.*(A+s.K-2)))));
%
% static autocorrelation term (norm)
static_term_norm = 2.*exp(-s.K).*(1-exp(s.K).*(1-s.K)).*(1-s.r)./s.K.^2;
%