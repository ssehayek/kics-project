% kicsSim(...) simulates two fluoroscent particle populations in 2D: one
% immobile population, which can be arranged along filaments or within a
% user-specified mask, and one freely diffusing population. Both
% populations are assumed to have the same photophysical properties (i.e.,
% blinking and bleaching rates). EMCCD noise is added to the simulations.
%
% The filaments are created by drawing angles from a narrow Gaussian
% distribution. Particles are then subsequently placed along the filaments
% with some probability. Distances between these particles are stochastic.
% See "directedFilaments.m" for more information.
%
% A fraction of the immobile population will be aggregated. The number of
% molecules in a certain aggregate is drawn from a Poisson distribution,
% with mean mean_agg_num. The distances between molecules in an aggregate
% and its center are drawn from a Gaussian distribution with standard
% deviation std_agg_dist.
%
% Please refer to "simulation-wrapper/kicsSimParams.m" for more details
% about the parameters.
%
% SIMULATION PARAMS
% sz: size of each frame in pixels (assumed square)
% T: number of frames
% w0: PSF size (e-2 radius)
% N_diff: number of diffusing molecules
% D: diffusion coefficient
% k_on: rate of blinking on (off -> on) in frames^-1
% k_off: rate of blinking off (on -> off) in frames^-1
% k_p: rate of bleaching in frames^-1
%
% AGGREGATE PARAMS
% prob_agg: probability that a certain immobile particle is aggregated
% mean_agg_num: mean number of particles in an aggregate
% std_agg_dist: std dev of distance between particles (from the initially
%               placed particle) - drawn from Gaussian of mean centered
%               around initially placed particle
%
% FILAMENT PARAMS
% num_filaments: number of filaments
% prob_place: probability to place molecule along filament - algorithm
%             steps along filament as: rand()*unit_vec, where unit_vec is
%             the unit vector in the direction of the filament
%
function [J,true_params] = kicsSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
    mean_agg_num,std_agg_dist,num_filaments,prob_place,varargin)

% assign defaults to varargin variables
%
frac_diff = [];
doSimulatePhotophysics = 1;
blink_model = 'twoStateBleach';
off_int_frac = 0;
n_sub_frames = 1;
save_run = 0;
use_parallel = 0;
n_parts = 1;
noise_type = '';
use_mask = 0;
kernel_varargin = {};
laser_varargin = {};
noise_varargin = {};
% assign user-specified varargin values to variables
for n = 1:2:length(varargin)
    if any(strcmpi(varargin{n},{'fracDiff','fractionDiff',...
            'fractionDiffusing'}))
        if isnumeric(varargin{n+1}) == 1
            frac_diff = varargin{n+1};
        else
            warning(['invalid option for varargin: ',varargin{n}]);
        end
    elseif any(strcmpi(varargin{n},{'parallel','useParallel'}))
        if length(varargin{n+1}) == 1 && any(varargin{n+1} == [0,1])
            use_parallel = varargin{n+1};
        else
            warning(['invalid option for varargin: ',varargin{n}]);
        end
    elseif any(strcmpi(varargin{n},{'doSimulatePhotophysics',...
            'simulatePhotophysics','photoSimulate','photoSim',...
            'simPhoto'})) && (varargin{n+1} == 0 || varargin{n+1} == 1)
        doSimulatePhotophysics = varargin{n+1};
    elseif any(strcmpi(varargin{n},{'model','blinkModel'}))
        if any(strcmpi(varargin{n+1},{'twoStateBleach','equalBleach',...
                'eqBleach'}))
        elseif any(strcmpi(varargin{n+1},{'offStateBleach','offBleach'}))
            blink_model = 'offStateBleach';
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'offIntFrac','partOffState'}))
        if isnumeric(varargin{n+1}) && 0 <= varargin{n+1} <= 1
            off_int_frac = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'nSubFrames','subFrames','subTimeSteps'})) && isnumeric(varargin{n+1})
        n_sub_frames = varargin{n+1};
    elseif any(strcmpi(varargin{n},{'savePath'}))
        if ischar(varargin{n+1})
            save_run = 1;
            savepath = varargin{n+1};
        else
            warning(['invalid option for varargin: ',varargin{n}]);
        end
    elseif any(strcmpi(varargin{n},{'useParallel','parallel'}))
        if any(varargin{n+1}==[0,1])
            use_parallel = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'nParts','numParts'}))
        if isnumeric(varargin{n+1}) && floor(varargin{n+1}) == varargin{n+1}
            n_parts = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'noiseType','noise'}))
        if any(strcmpi(varargin{n+1},{'emccd'}))
            noise_type = 'emccd';
        elseif any(strcmpi(varargin{n+1},{'legacy'}))
            noise_type = 'legacy';
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'mask','useMask'}))
        if isstruct(varargin{n+1}) && numel(fieldnames(varargin{n+1})) == 3
            mask_struct = varargin{n+1};
            
            use_mask = mask_struct.use_mask;
            mask_filepath = mask_struct.mask_filepath;
            N_imm = mask_struct.N_imm;
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'kernelVarargin','kernelVar','kerVar'}))
        if iscell(varargin{n+1})
            kernel_varargin = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'laserVarargin','laserVar'}))
        if iscell(varargin{n+1})
            laser_varargin = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{n},{'noiseVarargin','noiseVar'}))
        if iscell(varargin{n+1})
            noise_varargin = varargin{n+1};
        else
            warning(['Unknown option for ''',varargin{n},...
                ''', using default options.'])
        end
    end
end

% setting up save directory/file
if save_run
    [save_dir,~,save_ext] = fileparts(savepath);
    % prepend current directory path if savepath is empty
    if isempty(save_dir)
        savepath = [cd,filesep,savepath];
    end
    % append default extension '.mat' to savepath if extension is empty
    if isempty(save_ext)
        savepath = [savepath,'.mat'];
    end
    % check if file under same name already exists, throw an error if it does
    if exist(savepath,'file')
        error('file already exists; please specify another save path')
    end
end

% total number of frames (including sub frames)
total_T = n_sub_frames*T;
% interval between sub frames
sub_time = 1/n_sub_frames;

%% initialize positions

disp('generating dye positions')

tic
%

% struct which stores positions of all particle states, as well as
% aggregation state
particles = struct('immobile',struct,'diffusing',struct);

% immobile particle distribution
if use_mask
    % mask distribution
    load(mask_filepath,'mask')
    agg_pos.position = generateMaskPositions(N_imm,mask,w0);
else
    % filament structure
    [filaments,agg_pos] = directedFilaments(sz,num_filaments,prob_place);
end
% generate aggregates on filament structure
[particles.immobile,imm_positions,N_imm] = generateAggPositions(agg_pos.position,...
    mean_agg_num,std_agg_dist);
% initial diffusing particle positions
if ~isempty(frac_diff)
    % if it isn't empty, use frac_diff to determine N_diff (overrides
    % any specified value for N_diff)
    N_diff = round(frac_diff./(1-frac_diff).*N_imm);
end
% use N_diff to determine number of diffusing particles
init_diff_positions = generatePositions(N_diff,sz,w0);
%
toc

%% photophysics

disp('simulating photophysics')

tic
%
% total number of particles
N = N_imm + N_diff;
if doSimulatePhotophysics
    % simulate photophysics for all dyes
    [photo_state,obs_state,tint_obs_state] = ...
        simulatePhotophysics(N,total_T,sub_time,k_on,k_off,k_p,'model',...
        blink_model,'offIntFrac',off_int_frac);
    % allocate photophysical traces to diffusing population
    particles.diffusing.photoState = photo_state(1:N_diff,:,:);
    particles.diffusing.obsState = obs_state(1:N_diff,:,:);
    %
    
    % allocate photophysical traces to immobile population
    %
    % cumulative sum over number of aggregates
    agg_csum = [0,cumsum(particles.immobile.nDyes)];
    % organize photophysical properties into particles struct;
    % particles.immobile.photoState{n} stores the photophysical states of
    % the particles in aggregate n
    for n = 1:length(agg_csum)-1
        particles.immobile.photoState{n} = ...
            photo_state(N_diff+1+agg_csum(n):N_diff+agg_csum(n+1),:,:);
        particles.immobile.obsState{n} = ...
            obs_state(N_diff+1+agg_csum(n):N_diff+agg_csum(n+1),:,:);
    end
    %
end
%
toc

%% image series creation (immobile)

disp('placing immobile particles')

tic
%
% array to store noiseless image series of immobile particles
J_imm = zeros(sz,sz,T);

% cell array for storing PSFs
kernel_info = cell(1,N_imm);
for n = 1:N_imm
    [kernel_info{n}.psf_kernel,kernel_info{n}.xcoor,kernel_info{n}.ycoor] = ...
        getImgKernel(J_imm,imm_positions(n,:),w0,kernel_varargin{:});
end
%

% time-integrated observed states of immobile particles
imm_tint_obs_state = tint_obs_state(:,:,N_diff+1:end);

% construct immobile image series
for t = 1:T
    % find indices of particles that are not off for the whole frame t
    on_inds = find(imm_tint_obs_state(t,1,:))';
    for n = on_inds
        % PSF coordinates of nth dye
        xcoor_n = kernel_info{n}.xcoor;
        ycoor_n = kernel_info{n}.ycoor;
        % corresponding PSF values
        psf_kernel_imm_n = kernel_info{n}.psf_kernel;
        % time-integrated image at time t
        J_imm(ycoor_n,xcoor_n,t) = J_imm(ycoor_n,xcoor_n,t) + ...
            imm_tint_obs_state(t,1,n).*sub_time.*psf_kernel_imm_n;
    end
end
%

clear imm_tint_obs_state

toc

%% image series creation (diffusing)

disp('placing diffusing particles')

tic
%
% array to store noiseless image series of diffusing particles
J_diff = zeros(sz,sz,T);

% split frames into subsets so that diff_positions can be saved and cleared
% periodically
[T_parts,n_parts] = partitionArr(T,n_parts);

if save_run
    % open new matfile for saving diffusing positions
    m = matfile([save_dir,filesep,'diff_positions.mat'],'Writable',true);
    % initialize diff_positions
    m.diff_positions = zeros(N_diff,2,T,n_sub_frames);
end

for p = 1:n_parts
    % split arrays into parts for periodic saving/clearing
    %
    % vector of frames in subset p
    T_vec_p = T_parts{p};
    % number of frames in subset p
    T_p = length(T_vec_p);
    
    % slice variables for parfor
    %
    J_diff_p = J_diff(:,:,T_vec_p);
    % observed photostates of diffusing dyes at each sub frame
    diff_obs_state_p = obs_state(1:N_diff,T_vec_p,:);
    % generate diff_positions with chosen n_parts
    [diff_positions_p,next_diff_position] = simulateDiffusion(N_diff,D,...
        init_diff_positions,T_p*n_sub_frames,sub_time);
    if n_sub_frames ~= 1
        %%% if using subframes
        if use_parallel
            % initialize array for parfor loop purposes
            J_temp = J_diff_p;
            parfor t = 1:T_p
                for n = 1:N_diff
                    % observed states of particle n at time t over all sub frames
                    diff_obs_state_n_t = diff_obs_state_p(n,t,:);
                    % get indices of sub frames where particle is emitting
                    on_inds = find(diff_obs_state_n_t ~= 0)';
                    diff_positions_n_t = squeeze(diff_positions_p(n,:,t,:));
                    for i_on = on_inds
                        diff_psf_kernel_temp = ...
                            getImgKernel(J_temp,diff_positions_n_t(:,i_on),...
                            w0,'legacy',use_parallel,kernel_varargin{:});
                        
                        % time-integrated image at time t
                        diff_kernel = diff_obs_state_n_t(i_on).*sub_time.*...
                            diff_psf_kernel_temp;
                        J_diff_p(:,:,t) = J_diff_p(:,:,t) + diff_kernel;
                    end
                end
            end
        else
            for t = 1:T_p
                for n = 1:N_diff
                    diff_obs_state_n_t = diff_obs_state_p(n,t,:);
                    on_inds = find(diff_obs_state_n_t ~= 0)';
                    diff_positions_n_t = squeeze(diff_positions_p(n,:,t,:));
                    for i_on = on_inds
                        [diff_psf_kernel_temp,xcoor_temp,ycoor_temp] = ...
                            getImgKernel(J_diff_p,diff_positions_n_t(:,i_on),...
                            w0,'legacy',use_parallel,kernel_varargin{:});
                        
                        % time-integrated image at time t
                        J_diff_p(ycoor_temp,xcoor_temp,t) = J_diff_p(ycoor_temp,xcoor_temp,t) + ...
                            diff_obs_state_n_t(i_on).*sub_time.*diff_psf_kernel_temp;
                    end
                end
            end
        end
    else
        %%% if not using subframes
        if use_parallel
            % initialize array for parfor loop purposes
            J_temp = J_diff_p;
            parfor t = 1:T_p
                for n = 1:N_diff
                    % observed state of particle n at time t
                    diff_obs_state_n_t = diff_obs_state_p(n,t);
                    if diff_obs_state_n_t
                        diff_positions_n_t = squeeze(diff_positions_p(n,:,t));
                        diff_psf_kernel_temp = ...
                            getImgKernel(J_temp,diff_positions_n_t,...
                            w0,'legacy',use_parallel,kernel_varargin{:});
                        
                        % time-integrated image at time t
                        diff_kernel = diff_obs_state_n_t.*sub_time.*...
                            diff_psf_kernel_temp;
                        J_diff_p(:,:,t) = J_diff_p(:,:,t) + diff_kernel;
                    end
                end
            end
        else
            for t = 1:T_p
                for n = 1:N_diff
                    % observed state of particle n at time t
                    diff_obs_state_n_t = diff_obs_state_p(n,t);
                    if diff_obs_state_n_t
                        diff_positions_n_t = squeeze(diff_positions_p(n,:,t));
                        [diff_psf_kernel_temp,xcoor_temp,ycoor_temp] = ...
                            getImgKernel(J_diff_p,diff_positions_n_t,...
                            w0,'legacy',use_parallel,kernel_varargin{:});
                        
                        % time-integrated image at time t
                        J_diff_p(ycoor_temp,xcoor_temp,t) = J_diff_p(ycoor_temp,xcoor_temp,t) + ...
                            diff_obs_state_n_t.*sub_time.*diff_psf_kernel_temp;
                    end
                end
            end
        end
    end
    % initial diff_positions to use for generating diff_positions in next
    % subset
    init_diff_positions = next_diff_position;
    J_diff(:,:,T_vec_p) = J_diff_p;
    if save_run
        % partial save of diff_positions
        m.diff_positions(:,:,T_vec_p,:) = diff_positions_p;
    end
    clear diff_positions_p diff_obs_state_p J_diff_p
    
    progBar(p,n_parts);
end
%
toc

% clear unneeded variables
clear photo_state obs_state J_temp diff_obs_state diff_positions

J_sig = J_imm + J_diff;

%% add noise

disp('generating noise')

tic
%
if strcmp(noise_type,'emccd')
    J = addEMCCDNoise(J_sig,laser_varargin,noise_varargin{:});
elseif strcmp(noise_type,'legacy')
    [J,shot_noise_ser,wg_noise_ser] = addNoise(J_sig,noise_varargin{:});
else
    % rename J_sig -> J
    J = J_sig;
    clear J_sig
end
%
toc

%% save simulation info

disp('saving simulation')

tic
%
% simulation parameters to compare with ACF fitted parameters
true_params = [D,k_on./(k_on+k_off),k_on+k_off,N_diff/N];

% convert image series, J, to a less memory intensive numeric type
if max(J(:)) <= intmax('int16') && strcmp(noise_type,'emccd')
    J = int16(J);
end

if save_run
    % save particles struct
    save([save_dir,filesep,'particles.mat'],'particles','-v7.3');
    clear particles
    
    % save all remaining variables to savepath
    save(savepath,'-v7.3')
    
    % remove memory intensive variables
    all_vars = who;
    exclude_vars = {'J_diff','J_imm','J_sig','tint_obs_state',...
        'obs_kernel_stat'};
    csvars = setdiff(all_vars,exclude_vars);
    
    % save 'compressed' simulation
    [cfpath,cfname,cfext] = fileparts(savepath);
    cfname = [cfname,'--compressed'];
    compressed_savepath = [cfpath,filesep,cfname,cfext];
    save(compressed_savepath,csvars{:});
end
%
toc