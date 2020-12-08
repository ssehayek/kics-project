function [positions_rearr,next_position] = ...
    simulateDiffusion(N,D,init_positions,total_frames,dt)

steps = sqrt(2*D*dt)*randn(N,2,total_frames-1);

positions = cumsum(cat(3,init_positions,steps),3);

if nargout == 2
    next_step = sqrt(2*D*dt)*randn(N,2,1);
    next_position = positions(:,:,end) + next_step;
end

frames = total_frames*dt;
% number of sub-frames per frame
sub_frames = 1/dt;

% rearrange state arrays into more convenient form
positions_rearr = zeros(N,2,frames,sub_frames);
for t = 1:frames
    positions_rearr(:,:,t,:) = positions(:,:,(t-1)*sub_frames+1:t*sub_frames);
end