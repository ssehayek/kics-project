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
% determines whether to subtract by the local temporal mean of the image
% series using a temporal window that starts at each respective time point
time_win = 0;
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
    elseif any(strcmpi(varargin{ii},{'timeWindow','timeWin','window'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} > 0 && ...
                length(varargin{ii+1})<=2
            % set windowing boolean to 1
            time_win = 1;
            % window size specification
            win_k = varargin{ii+1};
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

J = double(J);

size_y = size(J,1);
size_x = size(J,2); % note inverted order definition of x and y
T = size(J,3);

if use_time_fluct && ~time_win
    % subtract by temporal mean
    J_fluct = J-repmat(mean(J,3),[1,1,size(J,3)]);
elseif time_win
    % subtract by local temporal mean
    J_fluct = J-movmean(J,[0,win_k-1],3);
else
    % no fluctuations
    J_fluct = J;
end

if force_even
    J_mirror=[J_fluct,fliplr(J_fluct(:,1:end-1,:));...
        flipud(J_fluct(1:end-1,1:end-1,:)),rot90(J_fluct(1:end-1,:,:),2)];
    J_k = fft2(J_mirror);
else
    J_k = fft2(J_fluct); % Fourier transform image series (space)
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
    phi_k = r_k./r_k(:,:,norm_lag+1);
else % don't normalize
    phi_k = r_k;
end