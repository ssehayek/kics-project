% return guesses and bounds on ICS parameters before performing fit using
% "ICSFit.m".
% 
% GUESSES
%
% amp_guess: the amplitude is guessed by interpolating the peak of the raw
% averaged ICS function with a spline
%
% offset_guess: the offset is guessed by averaging the absolute value of
% all points which are not included in "sub_corr" input (with the (0,0) lag
% re-added)
% 
% w0_guess: the e-2 radius of the PSF is guessed by explicitly finding
% which point is closest to having an e-2 ratio with the interpolated
% amplitude (after subtracting both values by "offset_guess")
%
% BOUNDS
%
% amp_bound: the amplitude is bounded from above by the max value of all
% raw ICS amplitudes (notice this value is raised by the noise, making it a
% natural choice as an upper bound; the bound is made tighter when the
% noise is small, or zero)
%
% offset_bound: similar to method for "offset_guess", except with averaging
% absolute values
%
% the bound for the PSF is more heuristic and is chosen to be 10, which is
% an exaggerated value
%
function [params_guess,lb,ub] = getICSGuess(corr,sub_corr)

avg_sub_corr = mean(sub_corr,3);

[mid_y,mid_x] = getCtrPxl(corr);
[sub_mid_y,sub_mid_x] = getCtrPxl(sub_corr);

% "sub_corr_00" is "avg_sub_corr" with the (0,0) lag value re-added, in case it
% has been removed
sub_corr_00 = sub_corr;  
sub_corr_00(sub_mid_y,sub_mid_x,:) = corr(mid_y,mid_x,:);

size_corr = [size(corr,1),size(corr,2)];
area_corr = prod(size_corr);
area_sub_corr = numel(avg_sub_corr);

%% guesses

% remove (0,0) lag, in case it hasn't already been removed, and interpolate
% its value with all other points in "sub_corr" input using a spline
tbi_sub_corr_avg = avg_sub_corr; 
tbi_sub_corr_avg(sub_mid_y,sub_mid_x) = NaN; % to be interpolated
amp_guess = interp2(tbi_sub_corr_avg,sub_mid_y,sub_mid_x,'spline');

sum_corr = sum(sum(mean(corr,3)));
sum_subcorr_00 = sum(sum(mean(sub_corr_00,3)));
offset_guess = 1/(area_corr-area_sub_corr)*(sum_corr-sum_subcorr_00);

% 
w0_dist = abs((avg_sub_corr-offset_guess)/(amp_guess-offset_guess)-exp(-2));
[min_dist_y,min_dist_x] = find(w0_dist==min(min(w0_dist)),1,'first');
w0_guess = norm(abs([sub_mid_y,sub_mid_x]-[min_dist_y,min_dist_x])); 

% vector output of guesses, sorted using the convention in "ICSFit.m"
params_guess = [amp_guess,w0_guess,offset_guess];

%% bounds

amp_bound = max(max(mean(corr,3)));

sum_abs_corr = sum(sum(mean(abs(corr),3)));
sum_abs_subcorr_00 = sum(sum(mean(abs(sub_corr_00),3)));
offset_bound = 1/(area_corr-area_sub_corr)*(sum_abs_corr-sum_abs_subcorr_00);

% take double of "offset_bound" in bounds for assurance
lb = [eps,eps,-2*offset_bound];
ub = [2*amp_bound,10,2*offset_bound];
