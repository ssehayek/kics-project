function [photo_state,obs_state] = simulatePhotophysics(N,total_T,...
    dt,k_on,k_off,k_p,varargin)

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

for t = 2:total_T
    for i = 1:N
        r = rand();
        if photo_state(i,t-1) == 1 % if previous state is on
            c = cumsum(p(:,1));
            j = find(c > r,1,'first');
            if j == 1 % j is next state
                % on -> on
                photo_state(i,t) = 1;
                obs_state(i,t) = 1;
            elseif j == 2
                % on -> off
                photo_state(i,t) = 2;
                obs_state(i,t) = off_int_frac;
            elseif j == 3
                % on -> bleach
                % no change in observed state, since it is initialized as 0
                photo_state(i,t) = 3;
            end
        elseif photo_state(i,t-1) == 2 % if previous state is off
            c = cumsum(p(:,2));
            j = find(c > r,1,'first');
            if j == 1
                % off -> on
                photo_state(i,t) = 1;
                obs_state(i,t) = 1;
            elseif j == 2
                % off -> off
                photo_state(i,t) = 2;
                obs_state(i,t) = off_int_frac;
            elseif j == 3
                % off -> bleach
                photo_state(i,t) = 3;
            end
        end
    end
end

