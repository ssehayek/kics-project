function [photo_state,obs_state] = simulatePhotophysics(N,total_T,...
    sub_time,k_on,k_off,k_p)

K = k_on + k_off;

% transition probabilities for blinking + bleaching
p_1_1 = exp(-k_p*sub_time)*(k_on + exp(-K*sub_time)*k_off)/K; %
p_2_1 = exp(-k_p*sub_time)*(1-exp(-K*sub_time))*k_off/K; % P(2|1)
p_3_1 = 1 - exp(-k_p*sub_time);
p_1_2 = exp(-k_p*sub_time)*(1-exp(-K*sub_time))*k_on/K;
p_2_2 = exp(-k_p*sub_time)*(exp(-K*sub_time)*k_on + k_off)/K;
p_3_2 = 1 - exp(-k_p*sub_time);

p_1 = [p_1_1,p_2_1,p_3_1];
p_2 = [p_1_2,p_2_2,p_3_2];
%

% photo-state initialization
%
% photo-states: 1:on, 2:off, 3:bleached
%
photo_state = 3*ones(N,total_T);
photo_state(:,1) = 1 + binornd(1,k_off/K,[N,1]);
%
% observed state initialization
%
% observed states. 0:off, 1:on
%
obs_state = (photo_state == 1);
%

for t = 2:total_T
    for i = 1:N
        r = rand();
        if photo_state(i,t-1) == 1 % if previous state is on
            c = cumsum(p_1);
            j = find(c > r,1,'first');
            if j == 1 % j is next state
                photo_state(i,t) = 1;
                obs_state(i,t) = 1;
            elseif j == 2
                % no change in observed state, since it is initialized as 0
                photo_state(i,t) = 2; 
            elseif j == 3
                % no change in observed state, since it is initialized as 0
                photo_state(i,t) = 3; 
            end
        elseif photo_state(i,t-1) == 2 % if previous state is off
            c = cumsum(p_2);
            j = find(c > r,1,'first');
            if j == 1
                photo_state(i,t) = 1;
                obs_state(i,t) = 1;
            elseif j == 2
                photo_state(i,t) = 2;
            elseif j == 3
                photo_state(i,t) = 3;
            end
        end
    end
end
    
