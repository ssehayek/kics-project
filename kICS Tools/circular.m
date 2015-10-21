% Written by S.S.
%
% Code for circular averaging
%
% This code is meant to circularly average about the center low frequency
% component of the autocorrelation function. Frequencies which have
% equivalent frequency magnitude, |k|, will be averaged.
%
% r_k_norm: Fourier transformed function R(k,tau) which
% has been shifted to have low frequencies in center of array by fftshift.
%
% X: Defines half square frame side length which defines the averaging range. Do
% not exceed the dimension of r_k. (Should be updated for non-squares)
%
% [r_k_abs, kSqVector] = circular(r_k_norm,X,varargin) returns a kICS
% correlation function averaged over |k|^2. kSqVector is the corresponding
% |k|^2 vector. 

function [r_k_abs,kSqVector] = circular(r_k_norm,X,varargin)

% r_k_abs has in its first column the possible frequency magnitudes in
% ascending order (more than possible - should be worked on). The second
% column is meant to place all equivalent frequencies. The third column
% counts the number of equivalent magnitudes. unique is meant to remove
% frequency magnitudes that have been multi-counted.
%

pixel_size = 1;

for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'pixelSize','pixelSz','sizeOfPixel','pixelWidth'}))
        pixel_size = varargin{i+1};
    end
end

size_y = size(r_k_norm,1);
size_x = size(r_k_norm,2);
T = size(r_k_norm,3);

kSqVector = zeros((X+1)^2,1);

x = 0;

for i = 0:X
    for j = 0:X
        x = x+1;
        kSqVector(x) = (i/(size_x*pixel_size))^2+(j/(size_y*pixel_size))^2;
    end
end

kSqVector = unique(kSqVector,'rows','sorted');
kSqVectorSz = length(kSqVector);

r_k_abs = zeros(kSqVectorSz,T);
kSqCounts = zeros(kSqVectorSz,1); % number of times kSqVector(i) occurs in r_k_norm

% The spiral loop. Starts at the frequency k=(0,0) in r_k and works its 
% way out.
% For a 3x3 the spiral ordering is shown:

% 432 %%%
% 501 %03
% 678 %12

% Each entry will be appropriately placed into its respective frequency
% magnitude bin by the loop.

u = 0;
v = 0;
du = 0;
dv = -1;

i = 1;

if ~mod(size_x,2)
    ctr_pxl_x = size_x/2+1;
else
    ctr_pxl_x = (size_x+1)/2;
end    

if ~mod(size_y,2)
    ctr_pxl_y = size_y/2+1;
else
    ctr_pxl_y = (size_y+1)/2;
end    

while i <= (2*X)^2
    if (-X < u <= X) && (-X < v <= X);
        index = find(kSqVector == (u/(size_x*pixel_size))^2+(v/(size_y*pixel_size))^2);
        try
            r_k_abs(index,:) = r_k_abs(index,:)+transpose(squeeze(r_k_norm(ctr_pxl_y+v,ctr_pxl_x+u,:)));
        catch
            disp([])
        end
        kSqCounts(index) = kSqCounts(index) + 1;
    end
    if u == v || (u < 0 && u == -v) || (u > 0 && u == 1-v);
		[du,dv] = deal(-dv, du);
    end
    [u,v] = deal(u+du, v+dv);
    i = i+1;
end

kSqVector = (2*pi)^2.*kSqVector;

% The remove is meant to take away higher frequencies which are beyond the
% range of the defined square. (To be improved)

% size(r_k_abs,1)
% 
% remove = find(r_k_abs(:,3)==0)
% r_k_abs(remove,:) = [];
% 
% r_k_abs(1:10,3)

r_k_abs = r_k_abs./repmat(kSqCounts,1,T); % Dividing the sum of each |k|^2 value by the appropriate number of occurences