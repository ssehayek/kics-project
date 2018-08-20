% Created by H.B. on 11/12/2014
% This function takes in an image pixel array and adds noise (Ref. [1])
% (1) From light source photon emission
% (2) From the probability of creating a photoelectron
% (3) By the EM register due to dark current and clock-induced charge
% (4) By electron shifting from detector array to read-out array
% (5) From read-out due to signal amplification 

% The scaling of number of photons emitted per molecule has functional form
% given by: Eq. 24 of ref. [2].  
% numPhotons = (fluorophore brightness)*numMolecules
% (fluorophore brightness) = I*beta*eta_w*exposureTime
% I = laser excitation intensity
% beta = excitation probability (detector quantum yield?)
% eta_w = detection efficiency

% Alternatively, the number of photons detected by the camera is shown also 
% in reference [3] Eqs. 11-12
% Ball-park range, we get about 10^5 photons/second detected


% REFERENCES: 
% [1] Hirsch et al. (2013) A stochastic model for electron multiplication
% charge coupled devices - from theory to practice, PloS One 8(1):e53671
% [2] Chen (1999) The photon counting histogram in fluorescence fluctuation
% spectroscopy, Biophys. J, 77: 553-67
% [3] Wohland et al. (2001) The standard deviation in fluorescence 
% correlation spectroscopy, Biophys. J, 80:2987-99


% USAGE: 
% imageSeries: [x y t] dimensional image series array*
% gain: emccd camera gain 
% avgPhotons: mean number of arriving photons (lambda) to detector per sec
% stdIntenisty: calculated by sqrt(2 g^2 avgPhotons)
% adf: analogue to digital conversion factor 

% ASSUMPTIONS
% *Assuming input simulation imageSeries shows the spatial probability of
% photon arrival (e.g. drawn from a normal distribution)
% ***An assumption made here is that the probability of detecting zero or
% negative values of electrons is negligible (i.e. lambda >> 1). 


function noisyImageSeries = addEMCCDNoise(imageSeries,gain,avgPhotons,...
            quantumYield,clockInducedCharge,darkCurrent,exposuretime,adf,readOutNoise)
if(~exist('gain','var'))
    gain = 164; % set default gain 
end
if(~exist('avgPhotons','var'))
    avgPhotons = 1e5; % (default) avg # photons/ sec / molecule
end
if(~exist('sdtIntensity','var')) % this is used if we don't assume (***)
    readOutNoise = 54;%sqrt(2*gain^2*avgPhotons); % default standard error
end
if(~exist('adf','var'))
    adf = 12; % default analogue to digital conversion
end
if(~exist('quantumYield','var'))
    quantumYield = 0.9; % default analogue to digital conversion
end
if(~exist('clockInducedCharge','var'))
    clockInducedCharge = 12; % baseline minimum spurious charge # electrons
end
if(~exist('darkCurrent','var'))
    darkCurrent = 0.0023; % electrons / pixel / second (dark current)
end
if(~exist('exposuretime','var'))
    exposuretime = 2e-3; % default exposure time
end


% compute the average number of input electrons detected by emccd per pixel
lambda = (avgPhotons*quantumYield*imageSeries + darkCurrent)*exposuretime + clockInducedCharge; 

% degOfFreedom = 4;
% 
% % compute the approximate imageSeries with noise (assuming *** that the
% % probability of negative image counts is negligible)
% noisyImageSeries = gain/2/adf*ncx2rnd(degOfFreedom,lambda); % factor of 2 off?

% Poisson noise contribution to  noise (due to source and detector)
nie = round(poissrnd(lambda)); % number of input electrons to EM register 

% EM register contribution to noise (not assuming )
noe = gamrnd(nie,gain); % number of output electrons from EM register 

% Readout contribution to noise
nread = noe + readOutNoise*randn(size(noe)); 

% shift negative values to zero
nread(nread<=0) = 0; 

% Convert final number of "output electrons" to image counts + discretize
noisyImageSeries = floor(nread/adf); 


end