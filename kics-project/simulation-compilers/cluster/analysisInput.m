%% SYSTEM DIRECTORIES

% path to codes (recursive)
code_path = '/media/3TB-GPU-Drive/simon.sehayek/research-projects/code';

% save directory for simulations
base_dir = '/media/3TB-GPU-Drive/simon.sehayek/research-projects/kics-reboot';

tmp_dir = '/media/3TB-GPU-Drive/simon.sehayek/research-projects/kics-reboot/queued-jobs/2020-10-15-19h26m03s/';

%% SIMULATION PARAMETERS

% simulation file tag
sim_tag = '';

% main parameters
%
% size of image series in pixels (sz x sz)
sz = 128;
% number of frames
T = 2048;
% number of subframes per frame to emulate detector time-integration
n_sub_frames = 50;
% PSF e-2 radius (pixels)
w0 = 3;
% number of diffusing particles
N_diff = 500;
% fraction of diffusing particles to immobile ones; this overrides the
% value for N_diff, unless it is left empty
frac_diff = 0.5;
% diffusion coefficient
D = 0.01;
% blinking on-rate (frames^-1)
k_on = 1;
% blinking off-rate (frames^-1)
k_off = 0.7;
% bleaching rate (frames^-1)
k_p = 0;

% aggregate parameters
%
% mean number of aggregates per immobile particle
mean_agg_num = 2;
% standard deviation of position from aggregate center
std_agg_dist = 0.1;

%%% IMMOBILE PARTICLE DISTRIBUTION
%
% FILAMENTS
%
% number of filaments
num_filaments = 20;
% probability of placing particles along filaments (prob for particles to
% be placed every rand() along filament)
prob_place = 0.3;
%
% MASK
%
% logical for using mask option (overrides filament options); mask must be
% a logical array with dimensions sz x sz; particles are distributed
% uniformly over area where mask == 1
use_mask = 0;
% if using mask, specify path to mask file; mask file must contain a
% logical array called "mask"
mask_filepath = '/media/3TB-GPU-Drive/simon.sehayek/research-projects/kics-reboot/masks/dendritic_cell_128x128.mat';
% number of immobile particles to distribute in mask
N_imm = 5000;
%
%%%

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
