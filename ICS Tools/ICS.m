% Written by S.O.S.
%
% This function is meant to take an image series and calculate its normalized
% ICS autocorrelaation function.

function [autocorr,xi,eta] = ICS(imgser,varargin)

mean_type = 'spatial'; % type of mean to subtract and normalize by
use_WKT = 1; % use Wiener-Khinchin theorem (WKT)

for i = 1:length(varargin)
    if strcmpi(varargin{i},'meanType')
        if strcmpi(varargin{i+1},'spatial')
            mean_type = 'spatial';
        elseif strcmpi(varargin{i+1},'temporal')
            mean_type = 'temporal';
        else
            error(['Invalid mean type specified; specify either',...
                ' ''spatial'' or ''temporal''.'])
        end
    elseif strcmpi(varargin{i},'useWKT')
        if varargin{i+1} == 1 || varargin{i+1} == 0
            use_WKT = varargin{i+1};
        else
            error(['Invalid option for "use_WKT" specified; value should,'...
                ' be either 0, or 1.'])
        end
    end
end

size_y = size(imgser,1);
size_x = size(imgser,2); % note inverted order definition of x and y
T = size(imgser,3);

% compute mean image series from specified option
% use of repmat is convenient for subtracting, and normalizing
if strcmp(mean_type,'spatial')
    mean_imgser = repmat(mean(mean(imgser,1),2),[size_y,size_x,1]);
else
    mean_imgser = repmat(mean(imgser,3),[1,1,T]);
end
fluct_imgser = imgser - mean_imgser;

% make a grid corresponding to lag values
% treat even & odd cases separately
if ~mod(size_x,2)
    xlags = -size_x/2:(size_x/2-1);
else
    xlags = -(size_x-1)/2:(size_x-1)/2;
end
if ~mod(size_y,2)
    ylags = -size_y/2:(size_y/2-1);
else
    ylags = -(size_y-1)/2:(size_y-1)/2;
end
[xi,eta] = meshgrid(xlags,ylags);

if use_WKT
    % for efficiency of fft, pad with zeros to closest power of 2 in both
    % dimensions
    if floor(log2(size_y)) == log2(size_y) % check if "size_y" is power of 2
        size_pad_y = 2*size_y; % length of FFT in first dimension with padding
    else % otherwise find next power of 2 for padded length
        pow2 = nextpow2(2*size_y);
        size_pad_y = 2^pow2;
    end
    if floor(log2(size_x)) == log2(size_x)
        size_pad_x = 2*size_y;
    else
        pow2 = nextpow2(2*size_x);
        size_pad_x = 2^pow2;
    end
    
    fluct_k = fft2(fluct_imgser,size_pad_y,size_pad_x);
    
    raw_autocorr = ifft2(fluct_k.*conj(fluct_k),'symmetric'); % autocorrelation using WKT (not yet properly normalized)
    
    % center zero-lags
    raw_autocorr = fftshift(fftshift(raw_autocorr,1),2);
    
    % crop out zero-padding
    ctr_pad_y = size_pad_y/2+1; ctr_pad_x = size_pad_x/2+1;
    crop_y = ctr_pad_y + ylags; crop_x = ctr_pad_x + xlags;
    raw_autocorr = raw_autocorr(crop_y,crop_x,:);
else
    % brute force computation of autocorrelation, sometimes required for
    % spatially inhomogeneous data
    raw_autocorr = zeros(size_y,size_x,T);
    for y = ylags
        for x = xlags
            eta_i = y + abs(ylags(1)) + 1;
            xi_i = x + abs(xlags(1)) + 1;
            if (sign(x) >= 0 && sign(y) >= 0) || sign(x) == sign(y)
                raw_autocorr(eta_i,xi_i,:) = sum(dot(...
                    fluct_imgser(1:size_y-abs(y),1:size_x-abs(x),:),...
                    fluct_imgser(abs(y)+1:size_y,abs(x)+1:size_x,:),1),2);
            elseif sign(y) < 0 && sign(x) >= 0
                raw_autocorr(eta_i,xi_i,:) = sum(dot(...
                    fluct_imgser(abs(y)+1:size_y,1:size_x-abs(x),:),...
                    fluct_imgser(1:size_y-abs(y),abs(x)+1:size_x,:),1),2);
            else % sign(y) >= 0 && sign(x) < 0
                raw_autocorr(eta_i,xi_i,:) = sum(dot(...
                    fluct_imgser(1:size_y-abs(y),abs(x)+1:size_x,:),...
                    fluct_imgser(abs(y)+1:size_y,1:size_x-abs(x),:),1),2);
            end
        end
    end
end

% normalize by spatial mean (temporal mean would not make any sense since
% it is pixel dependent, but maybe spatio-temporal could be a good
% choice...)
sp_mean_imgser = repmat(mean(mean(imgser,1),2),[size_y,size_x,1]);
autocorr = repmat(1./((size_y-abs(eta)).*(size_x-abs(xi))),1,1,T).*...
    raw_autocorr./sp_mean_imgser.^2;