% directedFilaments: This function is meant to produce a matrix of ones
% representing filaments in the cytoskeleton. This function is written for
% a square array, but should be generalized. Filaments in this code have
% angles drawn from a Gaussian with a fairly narrow std dev (by default), so
% that filaments are fairly parallel.
%
% (OLD UPDATE: This is a modification of filaments.m as an angle is now stored on each
% pixel containing a filament in "filaments" array. molecule_pos [now struct "particles"] keeps track
% of all stationary, bound proteins on cytoskeleton.)
%
% NEW UPDATE: This is a modification of filaments2.m as filament angles and
% particle positions are now stored in the struct "particles". Angles are
% no longer stored in the "filaments" array.
%
% M: Number of filaments
% sz: size of input array (sz x sz)
% pr: particles are placed on the filaments probabilistically, "pr" is this
% correponding probability
%
% Issues: Square array input
%         There's an obvious discretization problem - rand() now multiplies
%           incrementing unit vector (fix?)
%         (Size of particles array)
function [filaments,particles,N_imm] = directedFilaments(sz,M,pr,varargin)
%

mean_theta = pi*rand(); % randomly generated mean angle of filament.
std_theta = pi/50; % standard deviation of Gaussian dist for drawing filament angles
for i = 1:length(varargin)
    if any(strcmpi(varargin{i},{'MeanTheta','MeanAngle','MeanFilamentAngle'})) && isnumeric(varargin{i+1}) && 0 < varargin{i+1} < pi
        mean_theta = varargin{i+1}; % predefine mean angle of filaments
    elseif any(strcmpi(varargin{i},{'StdTheta','StdAngle','StdFilamentAngle'})) && isnumeric(varargin{i+1})
        std_theta = varargin{i+1}; % override default std dev of angle dist
    end
end

filaments = zeros(sz,sz); % indicator array of filament complex

start_points = zeros(M,1,2); % This array will store all start points of the filaments in the first sheet and their edges in the second sheet. 1 left, 2 top, 3 right, 4 bottom
start_points(:,:,1) = 1+(sz-1)*rand([M,1]); % Randomly distribute points on borders
start_points(:,:,2) = unidrnd(4,[M,1]); % Allocate start points to edges picked at random

line_positions = zeros(M,2); % Temporarily stores all true positions of line
particles = struct('position',zeros(sz^2,2),'angle',zeros(sz^2,1)); % Struct to store particles positions, and underlying filament angles.
% Note the size can be greater than preallocated.

% choose correct angle according to edge picked. Index of theta_shift
% is equal to the edge number picked in the upcoming loop. Angle is defined
% from bottom edge (3), in ccw direction

theta_shift = zeros(1,4);

j=1; % Iterates over particle number
for i = 1:M % loop over filament index
    boundary = [1,0,sz,0;0,1,0,sz]; % Map 1:4 to boundary points. Columns represent boundary coordinate (with 0 for replacement)
    p = boundary(:,start_points(i,1,2)); % Picking the right boundary for each point
    p(p==0) = start_points(i,1,1); % Replacing undetermined coordinate (0) with staring position.
    filaments(round(p(1)),round(p(2))) = 1; % Boundary point determined and placed
    
    line_positions(i,1) = p(1); % Initial positions
    line_positions(i,2) = p(2);
    if rand()<pr
        particles.position(j,:) = [line_positions(i,1),line_positions(i,2)];
        particles.angle(j) = mean_theta;
        j = j+1;
    end
    
    loop = true; % so long as line doesn't surpass border, loop = true
    
    theta = mean_theta + std_theta*randn();
    
    theta_shift(1) = theta + pi; % top edge
    theta_shift(3) = theta; % bottom edge
    if theta > pi/2
        [theta_shift(2),theta_shift(4)] = deal(theta + pi,theta); % 2: left edge 4: right edge
    else
        [theta_shift(2),theta_shift(4)] = deal(theta,theta + pi);
    end
    
    edge = start_points(i,:,2);
    dl = [-sin(theta_shift(edge)),cos(theta_shift(edge))];
    
    while loop
        line_positions(i,:) = line_positions(i,:) + rand()*dl; % note x & y inversion (also down is +ve y)
        if round(line_positions(i,1))<1||round(line_positions(i,1))>sz||round(line_positions(i,2))<1||round(line_positions(i,2))>sz
            loop = false;
        else
            filaments(round(line_positions(i,1)),round(line_positions(i,2)),1) = 1;
            % filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
            if rand()<pr
                particles.position(j,:) = [line_positions(i,1) line_positions(i,2)];
                particles.angle(j) = mean_theta;
                j = j+1;
            end
        end
    end
end
particles.position = particles.position(any(particles.position,2),:); % Remove all zero rows
particles.angle = particles.angle(any(particles.angle,2));
N_imm = size(particles.position,1);

% figure() % show filaments
% colormap(pink)
% imagesc(filaments)
% 
% figure() % show particle positions
% image = zeros(sz);
% for i = 1:size(particles.position,1)
%     image(round(particles.position(i,1)),round(particles.position(i,2))) = 1;
% end
% imagesc(image)
% colormap(pink)