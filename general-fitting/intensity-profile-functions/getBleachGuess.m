% return guesses and bounds on ICS parameters before performing fit using
% "bleachFit.m".
% 
% GUESSES
%
% amp_guess: the amplitude is estimated using the first point of the
% intesity trace
%
% offset_guess: the offset is guessed by using the last point of the
% intensity trace
% 
% kp_guess: the bleaching rate is guessed by explicitly finding which point
% in the intesity trace is closest to having an e-1 ratio with the
% guessed amplitude (after subtracting both values by "offset_guess")
%
function [params_guess] = getBleachGuess(J,varargin)
% function [params_guess,lb,ub] = getBleachGuess(J,varargin)

fit_offset = 0;
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'fitOffset','offset'}))
        fit_offset = 1;
        if length(varargin) >= ii+1 && isnumeric(varargin{ii+1}) && ...
                any(varargin{ii+1}==[0,1])
            fit_offset = varargin{ii+1};
        end
    end
end

int_trace = mean(mean(J,1),2);

%% guesses

amp_guess = int_trace(1);

offset_guess = int_trace(end);

kp_dist = abs((int_trace-offset_guess)/(amp_guess-offset_guess)-exp(-1));
min_dist = find(kp_dist==min(kp_dist),1,'first');
kp_guess = 1/min_dist;

% vector output of guesses, sorted using the convention in "bleachFit.m"
if fit_offset
    params_guess = [amp_guess,kp_guess,offset_guess];
else
    params_guess = [amp_guess,kp_guess];
end

% %% bounds
% 
% if fit_offset
%     lb = [offset_guess,eps,eps];
%     ub = 2*[amp_guess,kp_guess,offset_guess];
% else
%     lb = [offset_guess,eps];
%     ub = 2*[amp_guess,kp_guess];
% end