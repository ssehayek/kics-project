% addNoise(...) gives the option to add read and shot noise.
%
% INPUT 
% J_sig: image series of expected number of photon counts from signal i.e. 
%        values range from [0,1]

function [J,shot_noise_ser,wg_noise_ser] = addNoise(J_sig,varargin)

read_noise = 0;
shot_noise = 0;
for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'meanNoise','addMeanNoise'}))
        % add read noise relative to the spatial-mean of the image series as a
        % function of time
        % use: varargin{i+1}=[mu_noise/mean_signal,snr]
        if length(varargin{i+1}) == 2 && isnumeric(varargin{i+1}) && all(varargin{i+1} >= 0)
            read_noise = 1;
            read_noise_type = 'mean';
            read_noise_params = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    elseif any(strcmpi(varargin{i},{'peakNoise','addPeakNoise'}))
        % add read noise relative to the highest intensity value of the image
        % series
        % use: varargin{i+1}=[mu_noise/peak_signal,psnr]
        if length(varargin{i+1}) == 2 && isnumeric(varargin{i+1}) && all(varargin{i+1} >= 0)
            read_noise = 1;
            read_noise_type = 'peak';
            read_noise_params = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    elseif any(strcmpi(varargin{i},{'noReadNoise'}))
        % no read noise
    elseif any(strcmpi(varargin{i},{'shotNoise','addShotNoise'}))
        % add shot noise
        % value is expected number of photons in the signal at PSF peak
        if isnumeric(varargin{i+1}) && varargin{i+1} >= 1
            shot_noise = 1;
            n_photons = varargin{i+1};
        else
            warning(['invalid option for varargin: ',varargin{i}]);
        end
    elseif any(strcmpi(varargin{i},{'noShotNoise','noPoissNoise'}))
        % no shot noise
    end
end

size_y = size(J_sig,1);
size_x = size(J_sig,2);

% add read noise
wg_noise_ser = 0; % read noise
if read_noise
    if strcmp(read_noise_type,'mean')
        [mu_noise_frac,snr] = deal(read_noise_params(1),read_noise_params(2));
        sp_mean = repmat(mean(mean(J_sig,1),2),[size_y,size_x]); % spatial mean of image series
        wg_noise_ser =  mu_noise_frac.*sp_mean + sp_mean/snr.*randn(size(J_sig)); % white Gaussian noise
    else
        [mu_noise_frac,psnr] = deal(read_noise_params(1),read_noise_params(2));
        peak_signal = max(J_sig(:));
        wg_noise_ser =  mu_noise_frac.*peak_signal + peak_signal/psnr*randn(size(J_sig)); % white Gaussian noise
    end
end

% add shot noise
shot_noise_ser = 0;
if shot_noise
    % image series of expected number of photon counts (signal)
    J_sig = J_sig.*n_photons;
    % image series with shot noise (signal)
    J_shot = poissrnd(J_sig);
    
    % pixel values are now integer photon counts
    wg_noise_ser = round(wg_noise_ser*n_photons);
    J_sig = round(J_sig);
    %
    
    shot_noise_ser = J_shot - J_sig;
end

J = J_sig + wg_noise_ser + shot_noise_ser;