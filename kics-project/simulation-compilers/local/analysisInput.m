%% RUN TYPE

% boolean for creating simulation with given parameters
createSim = 1;
% boolean for fitting the simulation
fitSim = 1;

%% SAVE DIRECTORIES

% save directory
base_dir = 'C:\Users\Simon\Dropbox (Wiseman Research)\research-projects\kics-project';

% boolean for saving simulation
save_sim = 1;
% simulation file name
save_simfile = [base_dir,filesep,'sim',filesep,'sim'];

% boolean for saving analysis
save_run = 0;
% analysis file name
save_runpath = [base_dir,filesep,'run'];
% directory for analysis saving
analysis_runpath = [save_runpath,filesep,'analysis'];

%% SIMULATION PARAMETERS

% main parameters
%
% size of image series in pixels (sz x sz)
sz = 64;
% number of frames
T = 2048;
% number of subframes per frame to emulate detector time-integration
n_sub_frames = 50;
% PSF e-2 radius (pixels)
w0 = 3;
% number of diffusing particles
N_diff = 500;
% diffusion coefficient
D = 1;
% blinking on-rate (frames^-1)
k_on = 0.5;
% blinking off-rate (frames^-1)
k_off = 0.4;
% bleaching rate (frames^-1)
k_p = 5e-4;

% aggregate parameters
%
% mean number of aggregates per immobile particle
mean_agg_num = 2;
% standard deviation of position from aggregate center
std_agg_dist = 0.1;

% filament parameters
%
% number of filaments
num_filaments = 20;
% probability of placing particles along filaments (prob for particles to
% be placed every rand() along filament)
prob_place = 0;

% enter "twoStateBleach" (default), or "offStateBleach"
blink_model = 'twoStateBleach';
% off-state intensity emission, relative to on-state (default 0)
off_int_frac = 0;

sim_parallel = 1;

% name/value pairs for getImgKernel.m varargin
%
% 'kernelSize'
% 'boundCond' | 'periodic', or 'none'
% 'kerType' | 'legacy', or 'integrate'
%
kernel_varargin = {'boundCond','periodic'};

% name/value pairs for addLaserProfile.m varargin (called in addEMCCDNoise.m;
% addNoise.m no longer supports laser profile integration)
%
% 'laserWidth' | value in [0,inf]
% 'laserShift' | real number
%
laser_varargin = {'laserWidth',2*sz};

% enter 'legacy', or 'emccd'; 'legacy' option with empty noise_varargin
% does not simulate noise in image series
noise_type = 'emccd';
% name/value pairs for 'legacy' noise
%
% 'meanNoise', 'peakNoise' | [mu_noise/signal,(p)snr]
% 'shotNoise' | n_photons
%
% varargin name options for 'emccd' noise (see addEMCCDNoise.m for more
% details and defaults)
%
% 'gain', 'avgPhotons', 'readNoise', 'adf', 'qYield', 'CIC',
% 'darkCurrent', 'intTime', 'autofluorPer'
noise_varargin = {'avgPhotons',5e3,'intTime',0.05,'autofluorPer',0.05};
% noise_varargin = {'peakNoise',[0,Inf]};

%% kICS COMPUTATION

% choose whether to Fourier interpolate while circularly averaging
do_interp = 0;

% number of angles to sample for fixed |k|
n_theta = 1000;

% bounds of |k|^2 for which kICS tau = 0 lag (unnormalized) is noise dominated
% (should appear flat).
% CHOOSE THESE VALUES CAREFULLY
ksq_min_noise = 10;
ksq_max_noise = 15;

%% kICS FIT

% parallel boolean
fit_parallel = 0;
% if parallel is true, choose number of start points to try
startPts = 1e5;
% point tolerance
tolX = 0;
% objective function tolerance
tolFun = 0;
% output manymins object
output_mins = 0;

% lags to fit (actual tau values)
tauVector = [1:5];
kSqMin = 0.01;
% put 'max' to fit entire range
kSqMax = 1.2;

% Format [D,k_on,k_off,diff_frac]
params_guess = [rand rand rand rand];
lb = eps*ones(1,4);
ub = [Inf 1 Inf 1];

%% kICS PLOT

% lags to plot (subset of tauVector)
plotTauLags = tauVector;
nPtsFitPlot = 1000;