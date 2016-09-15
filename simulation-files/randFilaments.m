% randFilaments: This function is meant to randomly generate filaments, on
% which particles are placed. This function returns the placed particle
% positions, as well as the angle of the filament on whihc they reside.
% Additionally, this function returns an array which shows the spatial
% distribution of the filaments. Simulation of fluorescent tracks of these
% particles should be written in a separate code. This function is written
% for a square array, but should be generalized.
%
% (OLD UPDATE: This is a modification of filaments.m as an angle is now
% stored on each pixel containing a filament in "filaments" array.
% molecule_pos [now struct "particles"] keeps track of all stationary,
% bound proteins on cytoskeleton.)
%
% NEW UPDATE: This is a modification of filaments2.m as filament angles and
% particle positions are now stored in the struct "particles". Angles are
% no longer stored in the "filaments" array.
%
% M: Number of filaments sz: size of input array (sz x sz) pr: particles
% are placed on the filaments probabilistically, "pr" is this correponding
% probability
%
% Issues: Square array input
%         There's an obvious discretization problem - rand() now multiplies
%           incrementing unit vector (fix?)
%         (Size of particles array)
function [filaments,particles] = randFilaments(sz,M,pr)
%

filaments = zeros(sz,sz); % indicator array of filament complex

start_points = zeros(M,1,2); % This array will store all start points of the filaments in the first sheet and their edges in the second sheet. 1 left, 2 top, 3 right, 4 bottom
start_points(:,:,1) = 1+(sz-1)*rand([M,1]); % Randomly distribute points on borders
start_points(:,:,2) = unidrnd(4,[M,1]); % Allocate start points to edges picked at random

line_positions = zeros(M,2); % Temporarily stores all true positions of line
particles = struct('position',zeros(sz^2,2),'angle',zeros(sz^2,1)); % Struct to store particles positions, and underlying filament angles.
% Note the size can be greater than preallocated.

j=1; % Iterates over particle number
for i = 1:M % loop over filament index
    boundary = [1,0,sz,0;0,1,0,sz]; % Map 1:4 to boundary points. Columns represent boundary coordinate (with 0 for replacement)
    p = boundary(:,start_points(i,1,2)); % Picking the right boundary for each point
    p(p==0) = start_points(i,1,1); % Replacing undetermined coordinate (0) with staring position.
    filaments(round(p(1)),round(p(2))) = 1; % Boundary point determined and placed
    
    line_positions(i,1) = p(1); % Initial positions
    line_positions(i,2) = p(2);
    theta = pi*rand(); % randomly generated angle of filament
    if rand()<pr
        particles.position(j,:) = [line_positions(i,1),line_positions(i,2)];
        particles.angle(j) = theta; 
        j = j+1;
    end
    
    loop = true; % so long as line doesn't surpass border, loop = true
        
    % choose correct increment according to edge picked
    if start_points(i,1,2) == 1
        dl = [sin(theta),cos(theta)];
    elseif start_points(i,1,2) == 2
        dl = [cos(theta),sin(theta)];
    elseif start_points(i,1,2) == 3
        dl = [-sin(theta),cos(theta)];
    else
        dl = [cos(theta),-sin(theta)];
    end
    
    while loop
        line_positions(i,:) = line_positions(i,:) + rand()*dl; % note x & y inversion (also down is +ve y)
        if round(line_positions(i,1))<1||round(line_positions(i,1))>sz||round(line_positions(i,2))<1||round(line_positions(i,2))>sz
            loop = false;
        else
            filaments(round(line_positions(i,1)),round(line_positions(i,2)),1) = 1;
            % filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
            if rand()<pr
                particles.position(j,:) = [line_positions(i,1) line_positions(i,2)];
                particles.angle(j) = theta;
                j = j+1;
            end
        end
    end
end
particles.position = particles.position(any(particles.position,2),:); % Remove all zero rows
particles.angle = particles.angle(any(particles.angle,2));

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
