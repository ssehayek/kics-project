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
%   'simulatePhotoPhysics' | values: (default 1) 0
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
function [J,sim_info] = dronpaSim(sz,T,w0,N_diff,D,k_on,k_off,k_p,...
    prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,varargin)

% num_filaments = 20;
% prob_place = 0.3;
% std_agg_dist = 1/10; % default standard deviation of aggregates' distance from initially placed particle
doSimulatePhotophysics = 1;
addNoise = 0;
sub_frames = 1;
for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'doSimulatePhotophysics','simulatePhotophysics','photoSimulate','photoSim','simPhoto'})) && (varargin{i+1} == 0 || varargin{i+1} == 1)
        doSimulatePhotophysics = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'addNoise','noise','whiteNoise','snr'})) && isnumeric(varargin{i+1})  % add noise with varargin{i+1}=SNR 
        if isnumeric(varargin{i+1}) && varargin{i+1} > 0
            addNoise = 1;
            snr = varargin{i+1};    
        elseif any(strcmpi(varargin{i+1},{'noNoise','noiseless','none'}))
            addNoise = 0;
        else
            warning('invalid option for varargin "snr"');
        end
    elseif any(strcmpi(varargin{i},{'nSubFrames','subFrames','subTimeSteps'})) && isnumeric(varargin{i+1})
        sub_frames = varargin{i+1};
    end
end

% total number of frames (including sub frames)
total_T = sub_frames*T;
% interval between sub frames
sub_time = 1/sub_frames;

%% initialize filament structure

[~,immobile_pos] = directedFilaments(sz,num_filaments,prob_place);
N_imm = size(immobile_pos.position,1);

%% initialize particles array

particles = struct('immobile',struct,'diffusing',struct); % struct which stores positions of all particle states, as well as aggregation state
particles.immobile = struct('position',zeros(N_imm,2),'photoState',zeros(N_imm,total_T),...
    'obsState',zeros(N_imm,total_T),'aggState',zeros(N_imm,1),'aggregate',struct);
particles.immobile.aggregate = struct('position',{{}},'photoState',{{}},'obsState',{{}}); % aggregate positions and photo states
particles.diffusing = struct('position',zeros(N_diff,2,total_T),'photoState',zeros(N_diff,total_T),'obsState',zeros(N_diff,total_T));

particles.immobile.position = immobile_pos.position;
particles.immobile.aggregate.position = cell(N_imm,1);
particles.immobile.aggregate.photoState = cell(N_imm,1);
particles.immobile.aggregate.obsState = cell(N_imm,1);

particles.diffusing.position(:,:,1) = sz*rand(N_diff,2);

%% aggregated states

agg_state = binornd(1,prob_agg,[N_imm,1]); % random draws from binomial dist to determine aggregation state of each immobile particle
particles.immobile.aggState = agg_state; % 0:non-aggregated, 1:aggregated

agg_pos = [];
for n = 1:N_imm
    if particles.immobile.aggState(n)
        num_agg_n = poissrnd(mean_agg_num); % draw number of molecules in nth aggregate
        if num_agg_n ~= 0
            mean_particle_pos = particles.immobile.position(n,:); % for simplicity
            particles.immobile.aggregate.position{n} = repmat(mean_particle_pos,[num_agg_n,1]) + std_agg_dist*randn([num_agg_n,2]); % place aggregates with Gaussian dist
            
            add_pos = particles.immobile.aggregate.position{n};
            agg_pos = cat(1,agg_pos,add_pos); % concatenating all aggregate positions (not including initial aggregate)
        end
    end
end
N_agg = size(agg_pos,1); % number of particles in aggregated state (not including initial aggregate)
N_stat = N_imm + N_agg; % number of static particles (including aggregates)
N = N_imm + N_diff + N_agg; % total number of particles

%% photophysics

if doSimulatePhotophysics
    disp('simulating photophysics')
    tic
    
    % diffusing particles
    [particles.diffusing.photoState,particles.diffusing.obsState] = ...
        simulatePhotophysics(N_diff,total_T,sub_time,k_on,k_off,k_p);
    
    % immobile particles
    [particles.immobile.photoState,particles.immobile.obsState] = ...
        simulatePhotophysics(N_imm,total_T,sub_time,k_on,k_off,k_p);
    
    % aggregates
    [agg_photo_state,agg_obs_state] = ...
        simulatePhotophysics(N_agg,total_T,sub_time,k_on,k_off,k_p);
    cum_agg = 0; % cumulative number of aggregates
    for n = 1:N_imm
        if particles.immobile.aggState(n) ~= 0
            % number of aggregates in nth immobile particle
            num_agg_n = size(particles.immobile.aggregate.position{n},1); 
            if num_agg_n ~= 0
                particles.immobile.aggregate.photoState{n} = ...
                    agg_photo_state(cum_agg+1:cum_agg+num_agg_n,:);
                particles.immobile.aggregate.obsState{n} = ...
                    agg_obs_state(cum_agg+1:cum_agg+num_agg_n,:);
                
                cum_agg = cum_agg + num_agg_n;
            end
        end
    end
    toc
    %
end

%% position update
disp('simulating diffusion')
tic

for t = 2:total_T
    for i = 1:N_diff
        % propagate positions by diffusion
        particles.diffusing.position(i,:,t) = particles.diffusing.position(i,:,t-1) + sqrt(2*D*sub_time)*randn(1,2); 
    end
end

toc

%% Image series creation

% concatenate immobile molecules' positions
imm_positions = cat(1,particles.immobile.position,agg_pos);
imm_observed_state = cat(1,particles.immobile.obsState,agg_obs_state);

J = zeros(sz,sz,T);

%
disp('generating image series')
tic

%
disp('placing immobile molecules')

img_kernel_stat = zeros(sz,sz,N_stat);
for i = 1:N_stat
    % get ith immobile particle position
    pos = [imm_positions(i,1),imm_positions(i,2)];
    % corresponding img kernel 
    img_kernel_stat(:,:,i) = getImgKernel(J,pos,w0);
end

frame = 0;
for t = 1:total_T
    if mod(t-1,sub_frames) == 0
        frame = frame + 1;
    end
    
    J_t_stat = zeros(sz);
    for i = 1:N_stat
        if imm_observed_state(i,t) ~= 0
            % add particle signal to sub frame
            J_t_stat = J_t_stat + sub_time*imm_observed_state(i,t)*img_kernel_stat(:,:,i);            
        end
    end
    % add sub frame to whole frame
    J(:,:,frame) = J(:,:,frame) + J_t_stat;
end

toc
%

%
disp('placing diffusing molecules')
prog_bar = 1; % initialize progress status

frame = 0;
for t = 1:total_T
    if mod(t-1,sub_frames) == 0
        frame = frame + 1;
    end
    
    pos_x_t = particles.diffusing.position(:,1,t);
    pos_y_t = particles.diffusing.position(:,2,t);
    J_t = zeros(sz);
    for i = 1:N_diff
        if particles.diffusing.obsState(i,t) ~= 0
            % get ith diffusing particle position
            pos = [pos_y_t(i),pos_x_t(i)];
            % corresponding kernel
            img_kernel = getImgKernel(J,pos,w0);

            J_t = J_t + sub_time*particles.diffusing.obsState(i,t)*img_kernel;
        end
    end
    J(:,:,frame) = J(:,:,frame) + J_t;
    if frame/T >= prog_bar*0.1 % display updated progress every 10%
        disp(['progress: ',num2str(frame/T*100),'% after ',num2str(toc),' s.'])
        prog_bar = prog_bar + 1;
    end
end

toc
%

if addNoise
    for t = 1:T
        sp_mean_t = mean(mean(J(:,:,t),1),2); % spatial mean of image series at time t
        wgn_t = sp_mean_t/snr*randn(sz); % white Gaussian noise 
        J(:,:,t) = J(:,:,t) + wgn_t; % add noise to image at time t
    end
end

%% save simulation info
% store raw info in "sim_info"
sim_info.N_diff = N_diff; % number of diffusing particles (input param)
sim_info.N_imm = N_imm; % number of immobile particles
sim_info.N_agg = N_agg; % number of aggregated particles (disjoint from N_imm)
sim_info.N = N; % total number of particles
sim_info.mean_imgser = mean(mean(mean(J))); % mean over time and space of image series
sim_info.sub_frames = sub_frames;

% definitions for true parameters
I0 = 1; % PSF amplitude assumed to be 1 in some arbitrary units
b = 1; % brightness assumed to be 1 in some arbitrary units
V = sz^2; % volume
K = k_on + k_off; % sum of blinking rates
noise_factor = 4*V/(k_on^2/K^2*b^2*I0^2*w0^4*pi^2); % factor in front of noise term

% true fit params
sim_info.D = D; % diffusion
sim_info.k_on = k_on; % blink on rate
sim_info.k_off = k_off; % blink off rate
sim_info.k_p = k_p; % bleach rate
sim_info.f_d = N_diff/N; % diffusing fraction
sim_info.w0 = w0; % PSF size (e-2 radius)
if addNoise ~= 0
    sim_info.eta_p = noise_factor/N*(sim_info.mean_imgser/snr)^2; % noise term
else
    sim_info.eta_p = 0;
end

sim_info.true_params = [D,k_on,k_off,sim_info.f_d,w0,sim_info.eta_p]; % store all true fit params