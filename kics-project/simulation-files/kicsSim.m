% dronpaSim(...) is meant to simulate the behaviour of a certain Dronpa
% mutant bound to beta actin. The beta actin is allowed to freely diffuse,
% or remain fixed in filamentous structures. Immobile particles found on
% these filaments are also allowed to form aggregates. Additionally, all
% particles are assumed to possess the same blinking/bleaching
% distributions.
%
% The filaments are created by drawing angles from a narrow Gaussian
% distribution. Particles are then subsequently placed along the filaments
% with some probability. Distances between these particles are stochastic.
% See "directedFilaments.m" for more information.
%
% A fraction of the immobile population will be aggregated. The number of
% molecules in a certain aggregate is drawn from a Poisson distribution,
% with some fixed mean. The distances between molecules in an aggregate are
% drawn from a narrow Gaussian distribution.
%
% SIMULATION PARAMS
% sz: size of each frame (assumed square)
% T: number of frames
% w0: PSF size (e-2 radius)
% N_diff: number of diffusing molecules
% D: diffusion coefficient
% k_on: rate of blinking on
% k_off: rate of blinking off
% k_p: rate of bleaching
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
% VARARGIN
% - Choose whether to simulate photophysics
%   'doSimulatePhotophysics' | values: (default 1) 0
% - Add white-noise to simulation, with given signal-to-noise ratio (SNR)
%   (by default noise is not added)
%   'SNR' | value: value of (0,Inf), but more standard (1,~4]
% - Simulate integration time by choosing a number of sub frames for each
%   image. Ideally, the number of sub frames would be infinite.
%   'subFrames' | value: [1,Inf) (default 1 i.e. no integration time)
% FUTURE IMPROVEMENTS
% - Dronpa has odd blinking behaviour, which does not follow the standard
%   two-state model (see http://www.pnas.org/content/102/27/9511.full).
% - PSF should have stochastic size
% - Laser power should affect blinking on rate
%
function [J,true_params] = kicsSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
    mean_agg_num,std_agg_dist,num_filaments,prob_place,varargin)

% num_filaments = 20;
% prob_place = 0.3;
% std_agg_dist = 1/10; % default standard deviation of aggregates' distance from initially placed particle
doSimulatePhotophysics = 1;
blink_model = 'twoStateBleach';
off_int_frac = 0;
n_sub_frames = 1;
save_run = 0;
use_parallel = 0;
kernel_varargin = {};
laser_varargin = {};
noise_varargin = {};
for n = 1:2:length(varargin)
    if any(strcmpi(varargin{n},{'parallel','useParallel'}))
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
    elseif any(strcmpi(varargin{n},{'noiseType','noise'}))
        if any(strcmpi(varargin{n+1},{'emccd'}))
            noise_type = 'emccd';
        elseif any(strcmpi(varargin{n+1},{'legacy'}))
            noise_type = 'legacy';
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

% filament structure
[filaments,agg_pos] = directedFilaments(sz,num_filaments,prob_place);
% initial diffusing particle positions
particles.diffusing.position(:,:,1) = generatePositions(N_diff,sz,w0);
% generate aggregates on filament structure
[particles.immobile,imm_positions,N_imm] = generateAggPositions(agg_pos.position,...
    mean_agg_num,std_agg_dist);
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
    % organize photophysical properties into particles struct; organized
    % by aggregate
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

%% diffusion position update

disp('simulating diffusion')

tic
particles.diffusing.position = simulateDiffusion(N_diff,D,...
    particles.diffusing.position(:,:,1),total_T,sub_time);
diff_positions = particles.diffusing.position;
%
toc

%% image series creation (immobile)

disp('placing immobile particles')

tic
%
% array to store noiseless image series of immobile particles
J_imm = zeros(sz,sz,T);

% get dye PSF arrays and corresponding positions
if use_parallel
    % array for storing PSFs
    psf_kernel_imm = zeros(sz,sz,N_imm);
    for n = 1:N_imm
        psf_kernel_imm(:,:,n) = ...
            getImgKernel(J_imm,imm_positions(n,:),w0,'legacy',...
            use_parallel,kernel_varargin{:});
    end
else
    kernel_info = cell(1,N_imm);
    for n = 1:N_imm
        [kernel_info{n}.psf_kernel,kernel_info{n}.xcoor,kernel_info{n}.ycoor] = ...
            getImgKernel(J_imm,imm_positions(n,:),w0,'legacy',...
            use_parallel,kernel_varargin{:});
    end
end
%

% time-integrated observed states of immobile particles
imm_tint_obs_state = tint_obs_state(:,:,N_diff+1:end);

% construct immobile image series
if use_parallel
    for t = 1:T
        % find indices of particles that are not off for the whole frame t
        on_inds = find(imm_tint_obs_state(t,1,:));
        % multiply emitting dyes by their respective time-integrated
        % photostate
        obs_kernel_stat = imm_tint_obs_state(t,1,on_inds).*sub_time.*...
            psf_kernel_imm(:,:,on_inds);
        % time-integrated image at time t
        J_imm(:,:,t) = sum(obs_kernel_stat,3);
    end
    delete(gcp)
else
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
end
%
toc

%% image series creation (diffusing)

disp('placing diffusing particles')

tic
%
% array to store noiseless image series of diffusing particles
J_diff = zeros(sz,sz,T);

% observed photostates of diffusing dyes at each sub frame
diff_obs_state = obs_state(1:N_diff,:,:);

if use_parallel
    % initialize array for parfor loop purposes
    J_temp = J_diff;
    parfor t = 1:T
        for n = 1:N_diff
            % observed states of particle n at time t over all sub frames
            diff_obs_state_n_t = diff_obs_state(n,t,:);
            % isolate 
            on_inds = find(diff_obs_state_n_t ~= 0)';
            diff_positions_n_t = squeeze(diff_positions(n,:,t,:));
            for i_on = on_inds
                diff_psf_kernel_temp = ...
                    getImgKernel(J_temp,diff_positions_n_t(:,i_on),...
                    w0,'legacy',use_parallel,kernel_varargin{:});
                
                % time-integrated image at time t
                diff_kernel = diff_obs_state_n_t(i_on).*sub_time.*...
                    diff_psf_kernel_temp;
                J_diff(:,:,t) = J_diff(:,:,t) + diff_kernel;
            end
        end
    end
else
    for t = 1:T
        for n = 1:N_diff
            diff_obs_state_n_t = diff_obs_state(n,t,:);
            on_inds = find(diff_obs_state_n_t ~= 0)';
            diff_positions_n_t = squeeze(diff_positions(n,:,t,:));
            for i_on = on_inds
                [diff_psf_kernel_temp,xcoor_temp,ycoor_temp] = ...
                    getImgKernel(J_diff,diff_positions_n_t(:,i_on),...
                    w0,'legacy',use_parallel,kernel_varargin{:});
                
                % time-integrated image at time t
                J_diff(ycoor_temp,xcoor_temp,t) = J_diff(ycoor_temp,xcoor_temp,t) + ...
                    diff_obs_state_n_t(i_on).*sub_time.*diff_psf_kernel_temp;
            end
        end
    end
end
%
toc

J_sig = J_imm + J_diff;

%%

disp('generating noise')

tic
%
if strcmp(noise_type,'emccd')
    J = addEMCCDNoise(J_sig,laser_varargin,noise_varargin{:});
else
    [J,shot_noise_ser,wg_noise_ser] = addNoise(J_sig,noise_varargin{:});
end
%
toc

%% save simulation info

disp('saving simulation')

tic
%
clear photo_state obs_state imm_tint_obs_state J_temp diff_obs_state ...
    diff_positions psf_kernel_imm
true_params = [D,k_on,k_off,N_diff/N];
% prevent double saving of noiseless movies
if all(J_sig(:) == J(:))
    clear J_sig
end
if save_run
    [save_dir,~,save_ext] = fileparts(savepath);
    % prepend current directory path if savepath empty
    if isempty(save_dir)
        savepath = [cd,filesep,savepath];
    end
    % append default extension '.mat' to savepath if empty
    if isempty(save_ext)
        savepath = [savepath,'.mat'];
    end
    %
    if exist(savepath,'file')
        warning('file already exists; not overwritten')
    else
        save(savepath,'-v7.3')
    end
end
%
toc