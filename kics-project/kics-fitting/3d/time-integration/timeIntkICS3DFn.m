function [diff_term,static_term,diff_term_norm,static_term_norm] = ...
    timeIntkICS3DFn(s,tauGrid,kSqGrid)

A = s.D.*kSqGrid;

B1 = (A + s.K)./(4.*s.D);
B2 = A./(4.*s.D);

% factor of z0^2/(4*sqrt(pi)) is omitted, as it appears in all terms
C1 = s.r.*(1-s.r)./(2.*s.D).*exp(B1.*s.z0.^2);
C2 = s.r.^2./(2.*s.D).*exp(B2.*s.z0.^2);

% diffusing autocorrelation factor
diff_term_1 = C1.*(1/16).*B1.^(-3/2).*s.D.^(-1).*(2.*B1.^(1/2).*exp(1).^((-1).*B1.*( ...
    s.z0.^2+4.*s.D.*(1+tauGrid))).*(exp(1).^(8.*B1.*s.D).*(s.z0.^2+4.*s.D.*((-1)+tauGrid)).^( ...
    1/2)+(-2).*exp(1).^(4.*B1.*s.D).*(s.z0.^2+4.*s.D.*tauGrid).^(1/2)+(s.z0.^2+4.*s.D.* ...
    (1+tauGrid)).^(1/2))+pi.^(1/2).*(((-1)+2.*B1.*(s.z0.^2+4.*s.D.*((-1)+tauGrid))).* ...
    erf(B1.^(1/2).*(s.z0.^2+4.*s.D.*((-1)+tauGrid)).^(1/2))+(-2).*((-1)+2.*B1.* ...
    s.z0.^2+8.*B1.*s.D.*tauGrid).*erf(B1.^(1/2).*(s.z0.^2+4.*s.D.*tauGrid).^(1/2))+((-1)+2.* ...
    B1.*(s.z0.^2+4.*s.D.*(1+tauGrid))).*erf(B1.^(1/2).*(s.z0.^2+4.*s.D.*(1+tauGrid)).^(1/2)) ...
    ));

diff_term_2 = C2.*(1/16).*B2.^(-3/2).*s.D.^(-1).*(2.*B2.^(1/2).*exp(1).^((-1).*B2.*( ...
    s.z0.^2+4.*s.D.*(1+tauGrid))).*(exp(1).^(8.*B2.*s.D).*(s.z0.^2+4.*s.D.*((-1)+tauGrid)).^( ...
    1/2)+(-2).*exp(1).^(4.*B2.*s.D).*(s.z0.^2+4.*s.D.*tauGrid).^(1/2)+(s.z0.^2+4.*s.D.* ...
    (1+tauGrid)).^(1/2))+pi.^(1/2).*(((-1)+2.*B2.*(s.z0.^2+4.*s.D.*((-1)+tauGrid))).* ...
    erf(B2.^(1/2).*(s.z0.^2+4.*s.D.*((-1)+tauGrid)).^(1/2))+(-2).*((-1)+2.*B2.* ...
    s.z0.^2+8.*B2.*s.D.*tauGrid).*erf(B2.^(1/2).*(s.z0.^2+4.*s.D.*tauGrid).^(1/2))+((-1)+2.* ...
    B2.*(s.z0.^2+4.*s.D.*(1+tauGrid))).*erf(B2.^(1/2).*(s.z0.^2+4.*s.D.*(1+tauGrid)).^(1/2)) ...
    ));

diff_term = diff_term_1 + diff_term_2;
%
% static autocorrelation term
static_term = (-2).*exp(1).^((-1).*s.K.*tauGrid).*s.K.^(-2).*((-1)+s.r).*s.r.*((-1)+cosh(s.K))./s.z0;
%
% diffusing autocorrelation term (norm)
diff_term_norm_1 = C1.*(1/8).*B1.^(-3/2).*s.D.^(-1).*(2.*B1.^(1/2).*exp(1).^((-1).*B1.*(4.*s.D+ ...
  s.z0.^2)).*((-1).*exp(1).^(4.*B1.*s.D).*(s.z0.^2).^(1/2)+(4.*s.D+s.z0.^2).^( ...
  1/2))+pi.^(1/2).*(8.*B1.*s.D.*erfc(B1.^(1/2).*s.z0)+((-1)+2.*B1.*s.z0.^2).* ...
  erfc(B1.^(1/2).*(s.z0.^2).^(1/2))+(1+(-2).*B1.*(4.*s.D+s.z0.^2)).*erfc( ...
  B1.^(1/2).*(4.*s.D+s.z0.^2).^(1/2))));

diff_term_norm_2 = C2.*(1/8).*B2.^(-3/2).*s.D.^(-1).*(2.*B2.^(1/2).*exp(1).^((-1).*B2.*(4.*s.D+ ...
  s.z0.^2)).*((-1).*exp(1).^(4.*B2.*s.D).*(s.z0.^2).^(1/2)+(4.*s.D+s.z0.^2).^( ...
  1/2))+pi.^(1/2).*(8.*B2.*s.D.*erfc(B2.^(1/2).*s.z0)+((-1)+2.*B2.*s.z0.^2).* ...
  erfc(B2.^(1/2).*(s.z0.^2).^(1/2))+(1+(-2).*B2.*(4.*s.D+s.z0.^2)).*erfc( ...
  B2.^(1/2).*(4.*s.D+s.z0.^2).^(1/2))));

diff_term_norm = diff_term_norm_1 + diff_term_norm_2;
%
% static autocorrelation term (norm)
static_term_norm = (-2).*s.K.^(-2).*((-1)+exp(1).^((-1).*s.K)+s.K).*((-1)+s.r).*s.r./s.z0;
%