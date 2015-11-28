% d: diffusion coefficient
% kSq: k-squared vector
% a: on-rate
% b: off-rate
% T: total number frames
% f: ratio of static molecules to dynamic
% tau: time-lag
%
function out = biasFluctFunction(d,a,b,f,T,kSq,tau,varargin)

useGPU = 0;
for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'useGPU','GPU'}))
        if isnumeric(varargin{i+1})
            useGPU = varargin{i+1};
        end
    end
end

K = a + b;

[k,t] = meshgrid(sqrt(kSq),0:T-1); % all possible time lags

% disp(useGPU)
if useGPU
%     disp('transfering arrays to GPU...')
    k = gpuArray(k);
    t = gpuArray(t);
end

G_t = a/K^2*(a+b*exp(-K*t)); % photophysical autocorrelation
R_k_t = f*exp(-d*k.^2.*t)+1-f;

tau1 = 0:tau;
tau2 = 1:T-tau-1;
tau3 = tau+1:T-1;
tau4 = 1:T-1;

R1 = (T-tau).*sum(G_t(tau1+1,:).*R_k_t(tau1+1,:),1)+...
    sum(repmat(T-tau-tau2,[length(kSq),1])'.*G_t(tau2+1,:).*R_k_t(tau2+1,:),1)+...
    sum(repmat(T-tau3,[length(kSq),1])'.*G_t(tau3+1,:).*R_k_t(tau3+1,:),1);
R3 = T.*G_t(1,:).*R_k_t(1,:)+...
    2*sum(repmat(T-tau4,[length(kSq),1])'.*G_t(tau4+1,:).*R_k_t(tau4+1,:),1);

out = G_t(tau+1,:).*R_k_t(tau+1,:) + 1/T^2*R3 - 2*R1/(T*(T-tau));
