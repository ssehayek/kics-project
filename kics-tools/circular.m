% Written by S.O.S.
%
% Code for circular averaging
%
% This code is meant to circularly average about the center low frequency
% component of the autocorrelation function. Spatial frequencies which have
% equivalent frequency magnitude, |k|, will be averaged.
%
% r_k_norm: Fourier transformed function R(k,tau) which
% has been shifted to have low frequencies in center of array by fftshift.
%
% [r_k_circ, ksq] = circular(r_k_norm,varargin) returns a kICS correlation
% function averaged over |k|^2. First dimension of "r_k_circ" has values of
% "r_k_norm" after being circularly averaged over |k|^2. "ksq" is the
% corresponding |k|^2 vector. Second dimension represents the averaged
% value at fixed time-lag, which is an integer between [0,T-1].
%
function [r_k_circ,ksq] = circular(r_k_norm,varargin)

% compute size of input array over every dimension
size_y = size(r_k_norm,1);
size_x = size(r_k_norm,2);
T = size(r_k_norm,3);
%

% get 1-d vector of corresponding lattice points
% treat even & odd cases separately
if mod(size_x,2) == 0
	xgv = -size_x/2:size_x/2-1;
else
	xgv = -(size_x-1)/2:(size_x-1)/2;
end

if mod(size_y,2) == 0
	ygv = -size_y/2:size_y/2-1; 
else
	ygv = -(size_y-1)/2:(size_y-1)/2;
end
%

[X,Y] = meshgrid(xgv,ygv); % create xy-lattices with zeros at the center position
lattice_sqrd = (X/size_x).^2 + (Y/size_y).^2; % norm squared of lattice
lattice_nums = unique(lattice_sqrd); % unique values occuring in lattice sorted
                                     % into vector
l_nums = length(lattice_nums); % number of numbers in "lattice_nums"        
                                     
r_k_circ = zeros(l_nums,T); % array to fill with circular average of r_k_norm
unique_inds = zeros(l_nums,1);
for i = 1:length(lattice_nums) 
    n = lattice_nums(i); 
    inds = find(lattice_sqrd == n); % indices in "lattice_sqrd" with value "n"
    [subs_y,subs_x] = ind2sub(size(lattice_sqrd),inds); % convert linear indices to subscripts
    unique_inds(i) = inds(1);
    
	l_inds = length(inds); % number of indices
    for j = 1:l_inds % loop over individual xy-coordinates
        y = subs_y(j); x = subs_x(j); % fix coordinate (x,y)
        
        % average over values in "r_k_norm" with same normed spatial subscript 
        r_xy = reshape(r_k_norm(y,x,:),[1,T]); % take out singleton dim to match size(r_k_circ(i,:))
        r_k_circ(i,:) = r_k_circ(i,:) + r_xy/l_inds; % add contribution of "r_xy" to circular mean
		%
    end
end

ksq = (2*pi)^2*lattice_nums; % unique |k|^2 values sorted into vector