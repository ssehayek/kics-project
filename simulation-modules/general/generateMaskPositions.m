% generate positions within a logical mask (square dimensions)
%
% INPUT VARIABLES
% 
% N: number of particles within mask desired
% mask: logical mask array
% w0: PSF e-2 radius
%
function [positions] = generateMaskPositions(N,mask,w0)

positions = zeros(N,2);
sz = size(mask,1);

% loop index
n = 0;
while n < N
    position_n = (1 - w0) + (2*w0+sz-1).*rand(1,2);
    % catch out of bounds errors
    try
        if mask(round(position_n(1)),round(position_n(2)))
            % add coordinates to positions array if within mask
            n = n + 1;
            positions(n,:) = position_n;
        end
    catch
        continue
    end
end