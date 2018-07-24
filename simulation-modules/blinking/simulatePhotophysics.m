function [photo_state_rearr,obs_state_rearr,tint_obs_state] = ...
    simulatePhotophysics(N,total_T,dt,k_on,k_off,k_p,varargin)

blink_model = 'twoStateBleach';
off_int_frac = 0;
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'model','blinkModel'}))
        if any(strcmpi(varargin{ii+1},{'twoStateBleach','equalBleach',...
                'eqBleach'}))
            % blinking + bleaching from both states (default)
        elseif any(strcmpi(varargin{ii+1},{'offStateBleach','offBleach'}))
            blink_model = 'offStateBleach';
            % blinking + bleaching from off-state
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'offIntFrac','partOffState'}))
        if isnumeric(varargin{ii+1}) && 0 <= varargin{ii+1} <= 1
            off_int_frac = varargin{ii+1};
        end
    else
        warning(['Unknown varargin input ''',varargin{ii},'''.'])
    end
end

% sum of kinetic rates 
K = k_on + k_off;
% initialize array to store transition probabilities
p = zeros(3,2);
switch blink_model
    case 'twoStateBleach'
        % transition probabilities for blinking + two-state bleaching
        % 1:on, 2:off, 3:bleached
        %
        p(1,1) = exp(-k_p*dt)*(k_on + exp(-K*dt)*k_off)/K; %
        p(2,1) = exp(-k_p*dt)*(1-exp(-K*dt))*k_off/K; % P(2|1) or P(1->2)
        p(3,1) = 1 - exp(-k_p*dt);
        p(1,2) = exp(-k_p*dt)*(1-exp(-K*dt))*k_on/K;
        p(2,2) = exp(-k_p*dt)*(exp(-K*dt)*k_on + k_off)/K;
        p(3,2) = 1 - exp(-k_p*dt);
    otherwise
        % transition probabilities for blinking + off-state bleaching
        %
        K_p = k_on + k_off + k_p;
        D = sqrt(-4*k_p*k_off+K_p^2);
        
        p(1,1) = 1/(2*D)*exp(-dt*(K_p+D)/2)*(-k_p-k_on+k_off+D+exp(D*dt)*...
            (k_p+k_on-k_off+D));
        p(2,1) = 1/D*exp(-dt*(K_p+D)/2)*(exp(dt*D)-1)*k_off;
        p(3,1) = 0;
        p(1,2) = p(2,1)*k_on/k_off;
        p(2,2) = 1/(2*D)*exp(-dt*(K_p+D)/2)*(k_p+k_on-k_off+D+exp(D*dt)*...
            (-k_p-k_on+k_off+D));
        p(3,2) = 1 - p(1,2) - p(2,2);
end

% photo-state initialization
%
photo_state = 3*ones(N,total_T);
photo_state(:,1) = 1 + binornd(1,k_off/K,[N,1]);
%
% observed state initialization
% 0:off, 1:on
%
obs_state = (photo_state == 1);
% convert logical to double
obs_state = double(obs_state);
% assign off-state "intensity" to off_int_frac
obs_state(obs_state == 0) = off_int_frac;
%

% random numbers for state switching
r = rand(N,total_T-1);
% cumulative probabilities given initial photo-state is 1
c_1 = cumsum(p(:,1));
% initial photo-state is 2
c_2 = cumsum(p(:,2));
for t = 2:total_T
    % photo-states at t-1
    prev_photo_state = photo_state(:,t-1);
    % indices for which initial photo-state is 1 at t-1
    prev_state_1_inds = find(prev_photo_state == 1);
    % indices for which initial photo-state is 2 at t-1
    prev_state_2_inds = find(prev_photo_state == 2);
    % corresponding random numbers at these indices which will determine
    % next state
    r_state_1 = r(prev_state_1_inds,t-1);
    r_state_2 = r(prev_state_2_inds,t-1);
    
    % indices for which next state will be 1
    next_state_1_inds = union(prev_state_1_inds(c_1(1) > r_state_1),...
        prev_state_2_inds(c_2(1) > r_state_2));
    % indices for which next state will be 2
    next_state_2_inds = union(prev_state_1_inds(c_1(1) < r_state_1 & c_1(2) > r_state_1),...
        prev_state_2_inds(c_2(1) < r_state_2 & c_2(2) > r_state_2));
    
    % photo-states at t (all particles initialized as bleached)
    next_photo_state = photo_state(:,t);
    % assign next photo-states to corresponding indices
    next_photo_state(next_state_1_inds) = 1;
    next_photo_state(next_state_2_inds) = 2;
    % 
    photo_state(:,t) = next_photo_state;
    
    % observed states at t (all particles initialized as bleached)
    next_obs_state = obs_state(:,t);
    % assign next observed states to corresponding indices
    next_obs_state(next_state_1_inds) = 1;
    next_obs_state(next_state_2_inds) = off_int_frac;
    %
    obs_state(:,t) = next_obs_state;
end

%% time-integrated observed states

%%% added 2018-07-16; incompatible with older codes
%%%
% number of frames after time-integration
T = total_T*dt;
% number of sub-frames per frame
sub_frames = 1/dt;
% rearrange state arrays into more convenient form
obs_state_rearr = zeros(N,T,sub_frames);
photo_state_rearr = zeros(N,T,sub_frames);
for t = 1:T
    obs_state_rearr(:,t,:) = obs_state(:,(t-1)*sub_frames+1:t*sub_frames);
    photo_state_rearr(:,t,:) = photo_state(:,(t-1)*sub_frames+1:t*sub_frames);
end

% time-integrated observed particle state
% useful for later image series construction
tint_obs_state = sum(obs_state_rearr,3);
% permute dimensions for later multiplication with img_kernel_stat
tint_obs_state = permute(tint_obs_state,[2,3,1]);
%
%%%