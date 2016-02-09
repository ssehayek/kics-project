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
%
% FUTURE IMPROVEMENTS
% - Dronpa has odd blinking behaviour, which does not follow the standard
%   two-state model (see http://www.pnas.org/content/102/27/9511.full). 
% - PSF should have stochastic size
% - Laser power should affect blinking on rate
%
function [J,sim_info] = dronpaSim(sz,T,w0,N_diff,diffusion,k_on,k_off,k_p,...
    prob_agg,mean_agg_num,std_agg_dist,num_filaments,prob_place,varargin)

% num_filaments = 20;
% prob_place = 0.3;
% std_agg_dist = 1/10; % default standard deviation of aggregates' distance from initially placed particle
simulatePhotophysics = 1;
addNoise = 0;
for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'simulatePhotophysics','photoSimulate','photoSim','simPhoto'})) && (varargin{i+1} == 0 || varargin{i+1} == 1)
        simulatePhotophysics = varargin{i+1};
    elseif any(strcmpi(varargin{i},{'addNoise','noise','whiteNoise','snr'})) && isnumeric(varargin{i+1})  % add noise with varargin{i+1}=SNR 
        if isnumeric(varargin{i+1}) && varargin{i+1} > 0
            addNoise = 1;
            snr = varargin{i+1};    
        elseif any(strcmpi(varargin{i+1},{'noNoise','noiseless','none'}))
            addNoise = 0;
        else 
            disp('invalid option for varargin "snr"'); exit
        end
    end
end

addpath(genpath('C:\Users\SimonS\Dropbox (Personal)\Research\PhD\SOFI-Project'));

%% initialize filament structure

[~,immobile_pos] = directedFilaments(sz,num_filaments,prob_place);
N_imm = size(immobile_pos.position,1);

%% initialize particles array

particles = struct('immobile',struct,'diffusing',struct); % struct which stores positions of all particle states, as well as aggregation state
particles.immobile = struct('position',zeros(N_imm,2),'photoState',zeros(N_imm,T),...
    'obsState',zeros(N_imm,T),'aggState',zeros(N_imm,1),'aggregate',struct);
particles.immobile.aggregate = struct('position',{{}},'photoState',{{}},'obsState',{{}}); % aggregate positions and photo states
particles.diffusing = struct('position',zeros(N_diff,2,T),'photoState',zeros(N_diff,T),'obsState',zeros(N_diff,T));

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
        agg_num = poissrnd(mean_agg_num); % draw number of molecules in nth aggregate
        if agg_num ~= 0
            mean_particle_pos = particles.immobile.position(n,:); % for simplicity
            particles.immobile.aggregate.position{n} = repmat(mean_particle_pos,[agg_num,1]) + std_agg_dist*randn([agg_num,2]); % place aggregates with Gaussian dist
            
            add_pos = particles.immobile.aggregate.position{n};
            agg_pos = cat(1,agg_pos,add_pos); % concatenating all aggregate positions (not including initial aggregate)
        end
    end
end
N_agg = size(agg_pos,1); % number of particles in aggregated state (not including initial aggregate)
N = N_imm + N_diff + N_agg; % total number of particles

%% photo-state initialization

if simulatePhotophysics
    % photo-states. 1:on, 2:off, 3:bleached
    % observed states. 0:off, 1:on
    
    K = k_on + k_off;
    
    % diffusing particles
    particles.diffusing.photoState = 3*ones(N_diff,T);
    particles.diffusing.photoState(:,1) = 1 + binornd(1,k_off/K,[N_diff,1]);
    particles.diffusing.obsState(find(particles.diffusing.photoState == 1),1) = 1;
    
    % immobile particles
    particles.immobile.photoState = 3*ones(N_imm,T);
    particles.immobile.photoState(:,1) = 1 + binornd(1,k_off/K,[N_imm,1]);
    particles.immobile.obsState(find(particles.immobile.photoState == 1),1) = 1;
    
    % aggregates
%     agg_photo_state = [];
%     agg_obs_state = [];
    for n = 1:N_imm
        if particles.immobile.aggState(n) ~= 0
            agg_num = size(particles.immobile.aggregate.position{n},1); % find number of aggregates for nth particle
            if agg_num ~= 0
                particles.immobile.aggregate.photoState{n} = 3*ones(agg_num,T);
                particles.immobile.aggregate.photoState{n}(:,1) = 1 + binornd(1,k_off/K,[agg_num,1]); % draw initial on/off states uniformly
                particles.immobile.aggregate.obsState{n} = zeros(agg_num,T);
                particles.immobile.aggregate.obsState{n}(find(particles.immobile.aggregate.photoState{n}(:,1) == 1),1) = 1; % set "on" particles to have observed state of 1
                
%                 add_photo_state = particles.immobile.aggregate.photoState{n}(:,1); % for simplicity
%                 agg_photo_state = cat(1,agg_photo_state,add_photo_state);
%                 
%                 add_obs_state = particles.immobile.aggregate.obsState{n}(:,1); % for simplicity
%                 agg_obs_state = cat(1,agg_obs_state,add_obs_state);
            end
        end
    end
end


%% photo-state update

if simulatePhotophysics
    K = k_on + k_off;
    
    % transition probabilities for blinking + bleaching
    p_1_1 = exp(-k_p)*(k_on + exp(-K)*k_off)/K; %
    p_2_1 = exp(-k_p)*(1-exp(-K))*k_off/K; % P(2|1)
    p_3_1 = 1 - exp(-k_p);
    p_1_2 = exp(-k_p)*(1-exp(-K))*k_on/K;
    p_2_2 = exp(-k_p)*(exp(-K)*k_on + k_off)/K;
    p_3_2 = 1 - exp(-k_p);
    
    p_1 = [p_1_1,p_2_1,p_3_1];
    p_2 = [p_1_2,p_2_2,p_3_2];
    
    tic
    disp('simulating photophysics')
    
    % diffusing particles
    [photo_state,observed_state] = deal(particles.diffusing.photoState,particles.diffusing.obsState); % for simplicity
    for t = 2:T
        for i = 1:N_diff
            r = rand();
            if photo_state(i,t-1) == 1 % if previous state is on
                c = cumsum(p_1);
                j = find(c > r,1,'first');
                if j == 1 % j is next state
                    photo_state(i,t) = 1;
                    observed_state(i,t) = 1;
                elseif j == 2
                    photo_state(i,t) = 2; % no change in observed state, since it is initialized as 0
                elseif j == 3
                    photo_state(i,t) = 3; % no change in observed state, since it is initialized as 0
                end
            elseif photo_state(i,t-1) == 2 % is previous state is off
                c = cumsum(p_2);
                j = find(c > r,1,'first');
                if j == 1
                    photo_state(i,t) = 1;
                    observed_state(i,t) = 1;
                elseif j == 2
                    photo_state(i,t) = 2;
                elseif j == 3
                    photo_state(i,t) = 3;
                end
            end
        end
    end
    [particles.diffusing.photoState,particles.diffusing.obsState] = deal(photo_state,observed_state);
    
    % immobile particles
    [photo_state,observed_state] = deal(particles.immobile.photoState,particles.immobile.obsState);
    for t = 2:T
        for i = 1:N_imm
            r = rand();
            if photo_state(i,t-1) == 1
                c = cumsum(p_1);
                j = find(c > r,1,'first');
                if j == 1
                    photo_state(i,t) = 1;
                    observed_state(i,t) = 1;
                elseif j == 2
                    photo_state(i,t) = 2;
                elseif j == 3
                    photo_state(i,t) = 3;
                end
            elseif photo_state(i,t-1) == 2
                c = cumsum(p_2);
                j = find(c > r,1,'first');
                if j == 1
                    photo_state(i,t) = 1;
                    observed_state(i,t) = 1;
                elseif j == 2
                    photo_state(i,t) = 2;
                elseif j == 3
                    photo_state(i,t) = 3;
                end
            end
        end
    end
    [particles.immobile.photoState,particles.immobile.obsState] = deal(photo_state,observed_state);
    clear('photo_state','observed_state')
    
    % aggregates
    [agg_photo_state,agg_observed_state] = deal(particles.immobile.aggregate.photoState,particles.immobile.aggregate.obsState);
    for i = 1:N_imm
        agg_num = size(particles.immobile.aggregate.position{i},1);
        if agg_num ~= 0 % loop over number of particles in the ith aggregate
            for t = 2:T
                for n = 1:agg_num
                    r = rand();
                    if agg_photo_state{i}(n,t-1) == 1
                        c = cumsum(p_1);
                        j = find(c > r,1,'first');
                        if j == 1
                            agg_photo_state{i}(n,t) = 1;
                            agg_observed_state{i}(n,t) = 1;
                        elseif j == 2
                            agg_photo_state{i}(n,t) = 2;
                        elseif j == 3
                            agg_photo_state{i}(n,t) = 3;
                        end
                    elseif agg_photo_state{i}(n,t-1) == 2
                        c = cumsum(p_2);
                        j = find(c > r,1,'first');
                        if j == 1
                            agg_photo_state{i}(n,t) = 1;
                            agg_observed_state{i}(n,t) = 1;
                        elseif j == 2
                            agg_photo_state{i}(n,t) = 2;
                        elseif j == 3
                            agg_photo_state{i}(n,t) = 3;
                        end
                    end
                end
            end
        end
    end
    [particles.immobile.aggregate.photoState,particles.immobile.aggregate.obsState] = deal(agg_photo_state,agg_observed_state);
    toc
end

obs_state_temp = [];
for n = 1:N_imm
    obs_state_temp = cat(1,obs_state_temp,agg_observed_state{n});
end
agg_observed_state = obs_state_temp;
clear('obs_photo_state');

%% position update
disp('simulating diffusion')
tic

for t = 2:T
    for i = 1:N_diff
        particles.diffusing.position(i,:,t) = particles.diffusing.position(i,:,t-1) + sqrt(2*diffusion)*randn(1,2); % propagate positions by diffusion
    end
end

toc

%% Image series creation
if N_agg ~= 0
    positions = cat(1,repmat(particles.immobile.position,[1,1,T]),particles.diffusing.position,repmat(agg_pos,[1,1,T])); % positions of all particles WHICH ARE NOT AGGREGATED
    observed_state = cat(1,particles.immobile.obsState,particles.diffusing.obsState,agg_observed_state); % observed states of all particles WHICH ARE NOT AGGREGATED
else
    positions = cat(1,repmat(particles.immobile.position,[1,1,T]),particles.diffusing.position); % positions of all particles WHICH ARE NOT AGGREGATED
    observed_state = cat(1,particles.immobile.obsState,particles.diffusing.obsState); % observed states of all particles WHICH ARE NOT AGGREGATED
end
    
J = zeros(sz,sz,T);

kernelSize = ceil(3*w0); % makes Gaussian kernel with radius 3*PSFsize
% kernelSize = round((kernelSize-1)/2); %
[x,y] = meshgrid(-kernelSize:kernelSize,-kernelSize:kernelSize);
factor = 1; % before factor = 1/(2*pi*w0);

disp('generating image series')
tic
for t = 1:T
    for i = 1:N
        if observed_state(i,t) ~= 0
            dx = -round(positions(i,1,t))+positions(i,1,t); % dx shift from pixel center
            dy = -round(positions(i,2,t))+positions(i,2,t); % dy shift from pixel center
            
            arg =  -2*((x-dx).*(x-dx) + (y-dy).*(y-dy))/(w0^2);
            
            % create PSF Gaussian kernel
            kernel = exp(arg).*factor;
            kernel(kernel<eps*max(kernel(:))) = 0;
            %     sumk = sum(kernel(:));
            %     if sumk ~= 0,
            %       kernel  = kernel/sumk;
            %     end
            nonZeroEl = find(kernel);
            % find kernel index on image
            xcoor = mod(x(nonZeroEl) + round(positions(i,1,t)),sz);
            ycoor = mod(y(nonZeroEl) + round(positions(i,2,t)),sz);
            % fix MATLAB index of 1
            xcoor(xcoor==0) = sz;
            ycoor(ycoor==0) = sz;
            % make image
            imgKernel = full(sparse(ycoor,xcoor,kernel(nonZeroEl),sz,sz));
            % add image to video
            J(:,:,t) = J(:,:,t) + observed_state(i,t)*imgKernel;
        end
    end
end
% resulting image with rounded positions to place in appropriate pixel.
% Note modulo is taken for periodic boundaries.
plot(squeeze(sum(sum(J))))
toc
if addNoise
    for t = 1:T
        sp_mean_t = mean(mean(J(:,:,t),1),2);
        wgn_t = sp_mean_t/snr*randn(sz);
        J(:,:,t) = J(:,:,t) + wgn_t;
    end
end

% save simulation info
sim_info.N_diff = N_diff;
sim_info.N_imm = N_imm;
sim_info.N_agg = N_agg;
sim_info.N = N;