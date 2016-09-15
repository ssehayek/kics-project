% This function either outputs the least squares error, or the fit
% function. The fit function is for blinking fluorescing molecules which
% are both diffusing and stationary, it also assumes the k-space
% fluctuations were autocorrelated (and not just the image intensities),
% and that the PSF is Gaussian. Furthermore, the sampling nature of the
% autocorrelation is considered in the calculation of the fit function. The
% function was written to be able to also deal with image series containing
% white noise. Note the ACFs are fit globally i.e. simultaneously over all
% tau range specified by "tauVector".

function out = kICSFitBiasFluct(params,kSq,tauVector,T,varargin)

% default param values which may be changed by varargin
errBool = 0; % determines whether function outputs error, or fit function
tauNorm = 0; % time lag to normalize by (actual lag value)
noNorm = 0; % whether or not the kICS AC is normalized
useGPU = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'error','err','residual','res'})) 
        errBool = 1;
        ydata = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'tauNorm','normByTau','normTau','normByLag','normLag'}))
        if isnumeric(varargin{i+1})
            tauNorm = varargin{i+1};
        elseif any(strcmpi(varargin{i+1},{'none','noNorm',''})) 
            noNorm = 1;
        else
            disp('undefined varargin input for "normByTau" in "kICSFitBiasFluct"...'); exit
        end
    elseif any(strcmpi(varargin{i},{'useGPU','GPU'}))
        if isnumeric(varargin{i+1})
            useGPU = varargin{i+1};
        else
            disp('undefined varargin input for "useGPU" in "kICSFitBiasFluct"...'); exit;
        end
    end
end

% The field "A" is for a fudge factor which should be placed in the case of
% no normalization
if isnumeric(tauNorm)
    s = struct('diffusion',params(1),'k_on',params(2),'k_off',params(3),'frac',params(4),'w0',params(5),'sigma',params(6));
elseif noNorm == 1 && all(tauVector)
    s = struct('diffusion',params(1),'k_on',params(2),'k_off',params(3),'frac',params(4),'w0',params(5),'A',params(6));
elseif noNorm == 1 && ~all(tauVector)
    s = struct('diffusion',params(1),'k_on',params(2),'k_off',params(3),'frac',params(4),'w0',params(5),'sigma',params(6),'A',params(7));
else 
    disp('error in input of "tauVector", or varargin'); return
end

if useGPU
    F = zeros(length(kSq),length(tauVector),'gpuArray'); % best fit function
    err = zeros(1,length(tauVector),'gpuArray'); % least squares error
else
    F = zeros(length(kSq),length(tauVector)); % best fit function
    err = zeros(1,length(tauVector)); % least squares error
end

for tauInd = 1:length(tauVector)
    tau = tauVector(tauInd);
    
    Ik = exp(-s.w0^2/8*kSq); % Fourier transform of PSF
    f = biasFluctFunction(s.diffusion,s.k_on,s.k_off,s.frac,T,kSq,tau,'useGPU',useGPU);
    
    if ~all(size(Ik) == size(f)) % weird bug here, for some reason this condition is not always false... 
%         disp('bugged code patched...')
        Ik = Ik';
    end
    
    if noNorm == 0
        f_norm = biasFluctFunction(s.diffusion,s.k_on,s.k_off,s.frac,T,kSq,tauNorm,'useGPU',useGPU);
%         f = gather(f); f_norm = gather(f_norm);
        F(:,tauInd) = Ik.^2.*f./(Ik.^2.*f_norm+s.sigma);
    elseif tau ~= 0 && noNorm == 1 % no normalization + tau is non-zero 
        F(:,tauInd) = s.A*Ik.^2.*f; 
    elseif tau == 0 && noNorm == 1 % no normalization + tau is zero 
        F(:,tauInd) = s.A*Ik.^2.*f + s.sigma;
    end
        
    if errBool
        err(:,tauInd) = norm(ydata(:,tauVector(tauInd)+1) - F(:,tauInd));
    end
end

if ~errBool
    out = gather(F); % output best fit function
else
    err = sum(err); % sum LS error over all tau
    out = gather(err); % output 
end