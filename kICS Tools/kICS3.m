% Written by S.S.
%
% This function is meant to take an image series and calculate its normalized 
% kICS autocorrelaation function.

function [r_k_norm] = kICS3(J,varargin)

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

normByLag = 0;
noNorm = 0;
spatialPadFactor = 1;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'normByLag','normalizeByLag','normalizeByNthLag','tauLagNorm'}))
        if isnumeric(varargin{i+1})
            normByLag = varargin{i+1};
        elseif any(strcmpi(varargin{i+1},{'none','noNorm'}))
            noNorm = 1;
        end
    elseif any(strcmpi(varargin{i},{'spatialPad','spatialPadding','spPad'}))
        if isnumeric(varargin{i+1})
            spatialPadFactor = varargin{i+1};
        end
    end
end

size_y = spatialPadFactor*size(J,1);
size_x = spatialPadFactor*size(J,2); % note inverted order definition of x and y
T = size(J,3);
if floor(log2(T)) == log2(T)
    T_pad = 2*T;
else
    pow2 = nextpow2(2*T); 
    T_pad = 2^pow2;
end

r_k_norm = zeros(size_y,size_x,T);

J_k = fft2(J,size_y,size_x);

F_J_k = fft(J_k,T_pad,3);
r_k = ifft((F_J_k).*conj(F_J_k),[],3);

r_k(:,:,T+1:end) = [];

for tau = 0:T-1
    r_k(:,:,tau+1) = 1/(T-tau)*fftshift(r_k(:,:,tau+1));
    if ~noNorm
        r_k_norm(:,:,tau+1) = r_k(:,:,tau+1)./r_k(:,:,normByLag+1);
    else
        r_k_norm = r_k;
    end
end
