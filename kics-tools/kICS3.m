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
function [r_k_norm] = kICS3(J,varargin)

% time lag to normalize by
normByLag = 0;
% logical for normalization (0 will normalize)
noNorm = 0;
% use Wiener-Khinchin theorem; note this technically should
% not be used for non-stationary processes in time (e.g.
% photobleaching)
use_WKT = 1;
% extend images periodically in an even manner. This is recommended when
% the data is not intrinsically periodic across its boundaries (e.g. real
% data). Note using discrete cosine transform (DCT) is equivalent to this
% method, but is not used here to maintain convention and to enable user to
% perform FFT on result
force_even = 0;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'normByLag','normalizeByLag','normalizeByNthLag','tauLagNorm'}))
        if isnumeric(varargin{i+1})
            normByLag = varargin{i+1};
        elseif any(strcmpi(varargin{i+1},{'none','noNorm'}))
            noNorm = 1;
        end
    elseif any(strcmpi(varargin{i},{'useWKT','WKT'}))
        use_WKT = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'even','forceEven','mirrorMovie'}))
        force_even = 1;
    end
end

size_y = size(J,1);
size_x = size(J,2); % note inverted order definition of x and y
T = size(J,3);

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

r_k_norm = zeros(size(J_k)); % array to store kICS ACF
for tau = 0:T-1
    r_k(:,:,tau+1) = 1/(T-tau)*fftshift(r_k(:,:,tau+1));
    if ~noNorm
        r_k_norm(:,:,tau+1) = r_k(:,:,tau+1)./r_k(:,:,normByLag+1); % normalization
    else % don't normalize
        r_k_norm = r_k;
    end
end
