% simulate switching of two on/off states (1,2) with another on/off state
% (3) (topology shown below) with each state bleaching at the same rate 
% 1 <-> 3 <-> 2 
% each state can either be classified as "on", or "off"
%
% INPUT
%
% on_rates: [vector of length 2] kinetic rates into states 1 and 2,
% respectively
%
% off_rates: [vector of length 2] kinetic rates out of states 1 and 2,
% respectively
%
% state_type: [vector of length 3 | entries are 0, or 1] classifies states
% 1 to 3, respectively, as "on" states (1), or "off" states (0)
%
function [photo_state,obs_state] = simulateMultiPhotophysics(N,total_T,...
    dt,on_rates,off_rates,state_type,k_p,varargin)

off_int_frac = 0;
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'offIntFrac','partOffState'}))
        if isnumeric(varargin{ii+1}) && 0 <= varargin{ii+1} <= 1
            off_int_frac = varargin{ii+1};
        end
    end
end

% get on/off rates from input
k_on = on_rates(1);
k_on_2 = on_rates(2);
k_off = off_rates(1);
k_off_2 = off_rates(2);

% initial probabilities to be in state 1 and 2
p1 = k_on.*k_off_2./(k_on_2.*k_off+k_off_2.*(k_on+k_off));
p2 = k_on_2.*k_off./(k_on_2.*k_off+k_off_2.*(k_on+k_off));
p3 = 1-p1-p2;
% sum of blinking rates
K = k_on + k_off + k_on_2 + k_off_2;
% discriminant appearing in eigensystem
D = sqrt(k_on.^2+2*k_on.*(k_on_2+k_off-k_off_2)+(k_on_2-k_off+k_off_2).^2);

% relevant negative eigenvalues
L = -[-1/2*(K+D+2*k_p),1/2*(D-K-2*k_p),-k_p];

% columns of V are eigenvectors (relevant entries only)
V = [(K-2*k_off-D)./(2*(k_off-k_off_2)),-(K-2*k_off_2-D)./(2*(k_off-k_off_2)),1;...
    (K-2*k_off+D)./(2*(k_off-k_off_2)),-(K-2*k_off_2+D)./(2*(k_off-k_off_2)),1;...
    -p1,-p2,-(1-p1-p2)];
V = V';

% initial condition constants, C(i,j), where i denotes the initial condition
% i.e. for fixed i, C(i,j), or the ith row of C, is the set of constants for the solution with
% initial condition i
denom = 2.*D.*(k_on_2.*k_off+(k_on+k_off).*k_off_2);
C = zeros(3);
C(1,:) = [-k_off.*(k_on_2.*(2*k_off-k_off_2)+k_off_2.*(k_on+k_off-k_off_2+D)),...
    k_off.*(k_on_2.*(2*k_off-k_off_2)+k_off_2.*(k_on+k_off-k_off_2-D)),-denom]./denom;
C(2,:) = [-k_off_2.*(-k_on.*(k_off-2*k_off_2)+k_off.*(k_on_2-k_off+k_off_2+D)),...
    k_off_2.*(-k_on.*(k_off-2*k_off_2)+k_off.*(k_on_2-k_off+k_off_2-D)),-denom]./denom;
C(3,:) = [(k_on.^2.*k_off_2+k_on_2.*k_off.*(k_on_2-k_off+k_off_2+D)+...
    k_on.*(k_on_2.*(k_off+k_off_2)+k_off_2.*(k_off-k_off_2+D))),...
    (-k_on.^2.*k_off_2+k_on_2.*k_off.*(-k_on_2+k_off-k_off_2+D)+...
    k_on.*(-k_on_2.*(k_off+k_off_2)+k_off_2.*(-k_off+k_off_2+D))),...
    -denom]./denom;
% transition probabilities
%
% initialize array to store transition probabilities
p = zeros(4,3);
for i = 1:3
    for j = 1:3
        p(i,j) = sum(C(j,:).*V(i,:).*exp(-L.*dt),2); % P(i|j) or P(j->i)
    end
end
p(4,:) = 1 - sum(p,1);

% photo-state initialization
photo_state = 4*ones(N,total_T);
% draw multinomial random variables for initial state of system
% output is an array where each row is 0s and a 1 in the position of the
% drawn state
initial_states = mnrnd(1,[p1,p2,p3],N);
% convert each row to a number representing each state
initial_states = max(initial_states,[],2);
photo_state(:,1) = initial_states;
%
% observed state initialization
% 0:off, 1:on
%
obs_state = zeros(size(photo_state));
% loop through states to assign corresponding observed states
for s = 1:3
    obs_state(photo_state == s) = state_type(s)+(1-state_type(s))*off_int_frac;
end
%

for t = 2:total_T
    for i = 1:N
        r = rand();
        % execute if previous state of molecule is not bleached
        % otherwise, leave unaltered as all states are initialized as
        % bleached
        if photo_state(i,t-1) ~= 4 
            % find next state j
            c = cumsum(p(:,photo_state(i,t-1)));
            j = find(c > r,1,'first');
            photo_state(i,t) = j;            
            % map to corresponding observed state, j
            % bleached state is the default, so ignore this case
            if j ~= 4
                obs_state(i,t) = state_type(j)+(1-state_type(j))*off_int_frac;
            end
        end
    end
end
