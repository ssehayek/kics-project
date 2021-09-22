%% SYSTEM DIRECTORIES

% save directory for storing simulations
%
% character vector
base_dir = '';

% simulation file name
%
% 'simulation' (default) | character vector
sim_tag = '';

% logical for generating simulation in parallel
%
% 0 or 1
sim_parallel = 1;

% number of times to periodically save/clear diffusion positions (this
% variable is important for avoiding out of memory error)
%
% positive integer 
n_parts = 10;

%% SIMULATION PARAMETERS

% size of image series in pixels (sz x sz)
%
% positive integer
sz = 128;

% number of frames
%
% positive integer
T = 2048;

% number of subframes per frame to emulate detector time-integration
%
% positive integer
n_sub_frames = 50;

% PSF e-2 radius in pixels
%
% positive double
w0 = 3;

%%% Diffusion parameters %%%

% number of diffusing particles
%
% positive integer
N_diff = 500;

% fraction of diffusing particles (relative to total number of particles);
% this overrides the value for N_diff, unless it is left empty
%
% number between [0,1]
frac_diff = 0.5;

% diffusion coefficient
%
% positive double
D = 0.01;

%%% Photophysical parameters %%%

% blinking on-rate (off -> on) in frames^-1
%
% positive double
k_on = 1;

% blinking off-rate (on -> off) in frames^-1
%
% positive double
k_off = 0.7;

% bleaching rate in frames^-1
%
% positive double
k_p = 1e-4;

% bleaching model
%
% OPTIONS 
% "twoStateBleach": bleaching assumed to occur from both on/off fluorescence
% states with equal rate 
% "offStateBleach": bleaching assumed to occur from off fluorescence state
blink_model = 'twoStateBleach';

% off-state intensity emission, relative to on-state
%
% number between [0,1]
off_int_frac = 0;

%%% Aggregate parameters %%%

% mean number of aggregates per immobile particle (Poisson distributed)
%
% positive double
mean_agg_num = 2;

% standard deviation of position from aggregate center (normal distributed)
%
% positive double
std_agg_dist = 0.1;

%%% Filament parameters %%%

% number of filaments
%
% positive integer
num_filaments = 20;

% probability of placing particles along filaments (filaments are formed by
% incrementing by rand in the direction of the filament at each step; this
% is the probability that determines whether a particles is placed at each
% step)
%
% number between [0,1]
prob_place = 0.3;

%%% Mask parameeters (overrides filament options) %%%

% logical for using mask option; mask must be
% a logical array with dimensions sz x sz; particles are distributed
% uniformly over area where mask == 1
%
% 0 or 1
use_mask = 0;

% if using mask, specify path to mask file; mask file must contain a
% logical array called "mask"
%
% character vector
mask_filepath = '';

% number of immobile particles to distribute uniformly within mask
%
% positive integer
N_imm = 5000;

%%% Varargin (check files mentioned for more info) %%%

% name/value pairs for getImgKernel.m varargin
%
% cell array with the following name/value pairs:
% 'kernelSize' | positive double 
% 'boundCond' | 'periodic', or 'none'
% 'kerType' | 'legacy', or 'integrate'
kernel_varargin = {'boundCond','periodic'};

% name/value pairs for addLaserProfile.m varargin (called in addEMCCDNoise.m;
% addNoise.m no longer supports laser profile integration)
%
% cell array with the following name/value pairs:
% 'laserWidth' | positive double
% 'laserShift' | double
laser_varargin = {'laserWidth',2*sz};

% noise type
%
% OPTIONS
% 'legacy': adds shot noise and read noise (see addNoise.m)
% 'emccd'; adds EMCCD noise (see addEMCCDNoise.m)
noise_type = 'emccd';

% name/value pairs for 'legacy' noise (see addNoise.m)
%
% cell array with the following name/value pairs:
% 'meanNoise', 'peakNoise' | [mu_noise/signal,(p)snr]
% 'shotNoise' | n_photons
%
% varargin name options for 'emccd' noise (see addEMCCDNoise.m for more
% details and defaults)
%
% cell array with the following names:
% 'gain', 'avgPhotons', 'readNoise', 'adf', 'qYield', 'CIC',
% 'darkCurrent', 'intTime', 'autofluorPer'
noise_varargin = {'avgPhotons',5e3,'intTime',0.05,'autofluorPer',0.05};