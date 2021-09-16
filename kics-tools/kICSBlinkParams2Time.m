% get on/off-times and confidence intervals from extended kICS fit
%
% INPUT
%
% x: best-fit parameters with x(:,1)=rho_on and x(:,2)=K (in inverse
% frames)
%
% dt: time between frames
function [blink_times,blink_times_se] = kICSBlinkParams2Time(x,dt)

% number of measurements
n = size(x,1);

% compute blinking rates
k_on = x(:,1).*x(:,2);
k_off = (1-x(:,1)).*x(:,2);

% convert rates into times, blink_times=[t_on,t_off]
blink_times(1) = mean(1./k_off)*dt;
blink_times(2) = mean(1./k_on)*dt;

% standard error on blink_times
blink_times_se(1) = std(dt./k_off)/sqrt(n);
blink_times_se(2) = std(dt./k_on)/sqrt(n);