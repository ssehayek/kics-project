% Written by Simon Sehayek
%
% This function is meant to take an image series and calculate its normalized
% kICS autocorrelation function.
%
% Taking the image into Fourier space J_k is the result of applying the
% fft2 on the image. fft2 takes the Fourier Transform of each J(:,:,t). r_k
% is the normalized autocorrelation function. Once the time variable is
% "summed out" by the autocorrelation, we can finally apply fftshift to get
% the low frequencies at the center.
%
% The autocorrelation function is computed through the circular convolution
% of the sequence.
%
% Notice the zero-padding since there is a mismatch in the definition of
% autocorrelation and circular convolution (circular convolution theorem
% and cross-correlation theorem).
% https://www.mathworks.com/help/signal/ug/linear-and-circular-convolution.html.
%
function [phi_k] = kICS(J,varargin)

% time lag to normalize by
norm_lag = 0;
% if normalizing by 0th lag, this option first subtracts the noise from
% this lag past ksq_min_noise
sub_noise = 0;
% logical for normalization
use_norm = 1;
% use Wiener-Khinchin theorem; note this technically should
% not be used for non-stationary processes in time (e.g.
% photobleaching)
use_WKT = 1;
% determines whether to subtract by the temporal mean of the image series.
% Note that subtracting by the spatial mean would leave the spatial Fourier
% transform unchanged
use_time_fluct = 1;
% extend images periodically in an even manner. This is recommended when
% the data is not intrinsically periodic across its boundaries (e.g. real
% data). Note using discrete cosine transform (DCT) is equivalent to this
% method, but is not used here to maintain convention and to enable user to
% perform FFT on result
force_even = 0;
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'normByLag','normalizeByLag','normalizeByNthLag','tauLagNorm'}))
        if isnumeric(varargin{ii+1}) && any(varargin{ii+1} == [0,1])
            norm_lag = varargin{ii+1};
        elseif any(strcmpi(varargin{ii+1},{'none','noNorm'}))
            use_norm = 0;
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'ksqMinNoise'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} >= 0
            sub_noise = 1;
            ksq_min_noise = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'useWKT','WKT'}))
        if isnumeric(varargin{ii+1}) && any(varargin{ii+1} == [0,1])
            use_WKT = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'useTimeFluct','subTempMean','subMean'}))
        if isnumeric(varargin{ii+1}) && any(varargin{ii+1} == [0,1])
            use_time_fluct = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'even','forceEven','mirrorMovie'}))
        force_even = 1;
    else
        warning(['unknown varargin input ''',varargin{ii},'''.'])
    end
end

if ( norm_lag ~= 0 || use_norm == 0 ) && sub_noise == 1
    warning('Noise subtraction specified, but not normalizing by tau = 0. Noise will not be subtracted.')
end

size_y = size(J,1);
size_x = size(J,2); % note inverted order definition of x and y
T = size(J,3);

if use_time_fluct % subtract by temporal mean
    J = J-repmat(mean(J,3),[1,1,size(J,3)]);
end

if force_even
    J_mirror=[J,fliplr(J(:,1:end-1,:));flipud(J(1:end-1,1:end-1,:)),rot90(J(1:end-1,:,:),2)];
    J_k = fft2(J_mirror);
else
    J_k = fft2(J); % Fourier transform image series (space)
end

if use_WKT
    if floor(log2(T)) == log2(T) % check if "T" is power of 2
        T_pad = 2*T; % length of FFT in time with padding
    else
        pow2 = nextpow2(2*T); % otherwise find next power of 2 of padded length
        T_pad = 2^pow2;
    end
    F_J_k = fft(J_k,T_pad,3); % Fourier transform again (time)
    r_k = ifft((F_J_k).*conj(F_J_k),[],3); % autocorrelation from Wiener-Khinchin (time)
    
    r_k(:,:,T+1:end) = []; % remove temporal padding
else
    r_k = zeros(size_y,size_x,T);
    for tau = 0:T-1
        r_k(:,:,tau+1) = dot(J_k(:,:,1:T-tau),J_k(:,:,tau+1:T),3);
    end
end

% shift (kx,ky)=(0,0) lag to be in center
for tau = 0:T-1
    r_k(:,:,tau+1) = 1/(T-tau)*fftshift(r_k(:,:,tau+1));
end

if use_norm % normalize
    if sub_noise && norm_lag == 0
        [~,~,noise_inds] = getKSqVector(r_k,'ksqMin',ksq_min_noise);
        r_k_0 = r_k(:,:,1);
        mean_noise = mean(r_k_0(noise_inds));
        phi_k = r_k./(r_k_0-mean_noise);
    else
        phi_k = r_k./r_k(:,:,norm_lag+1);
    end
else % don't normalize
    phi_k = r_k;
end
