% randFilaments: This function is meant to produce a matrix of ones
% representing filaments in the cytoskeleton. This function is written for
% a square array, but should be generalized.
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
%         (Size of particles)
function [filaments,particles] = randFilaments(sz,M,pr)
%

filaments = zeros(sz,sz); % 1st sheet is actual image of cytoskeleton. 2nd sheet is the angle of filament.

start_points = zeros(M,1,2); % This array will store all start points of the filaments in the first sheet and their edges in the second sheet. 1 left, 2 top, 3 right, 4 bottom
start_points(:,:,1) = 1+(sz-1)*rand([M,1]); % Randomly distribute points on borders
start_points(:,:,2) = unidrnd(4,[M,1]); % Allocate start points to edges picked at random

line_positions = zeros(M,2,1); % Temporarily stores all true positions of line
particles = struct('position',zeros(sz^2,2),'angle',zeros(sz^2,1)); % Struct to store particles positions, and underlying filament angles. 
% Note the size can be greater than preallocated. 

j=1; % Iterates over particle number
for i = 1:M
    boundary = [1,0,sz,0;0,1,0,sz]; % Map 1:4 to boundary points. Columns represent boundary coordinate (with 0 for replacement)
    p = boundary(:,start_points(i,1,2)); % Picking the right boundary for each point
    p(p==0) = start_points(i,1,1); % Replacing undetermined coordinate (0) with staring position.
    filaments(round(p(1)),round(p(2))) = 1; % Boundary point determined and placed
    if rand()<pr 
        particles.position(j,:) = [line_positions(i,1,1),line_positions(i,2,1)];        
        j = j+1;
    end
    
    line_positions(i,1,1) = p(1); % Initial positions
    line_positions(i,2,1) = p(2);
    loop = true; % So long as line doesn't surpass border, loop = true
    
    % Here we treat the 4 different edge cases
    if start_points(i,1,2) == 1
        theta = pi*rand()
        while loop
            line_positions(i,1,1) = line_positions(i,1,1) + sin(theta); % note x & y inversion (also down is + y)
            line_positions(i,2,1) = line_positions(i,2,1) + cos(theta);
            if round(line_positions(i,1,1))<1||round(line_positions(i,1,1))>sz||round(line_positions(i,2,1))<1||round(line_positions(i,2,1))>sz
                loop = false;
            else
                filaments(round(line_positions(i,1,1)),round(line_positions(i,2,1)),1) = 1;
%                 filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
                if rand()<pr
                    particles.position(j,:) = [line_positions(i,1,1) line_positions(i,2,1)];
                    j = j+1;
                end
            end
        end
    elseif start_points(i,1,2) == 2
        theta = pi*rand();
        while loop
            line_positions(i,1,1) = line_positions(i,1,1) + cos(theta); % x & y inversion
            line_positions(i,2,1) = line_positions(i,2,1) + sin(theta);
            if round(line_positions(i,1,1))<1||round(line_positions(i,1,1))>sz||round(line_positions(i,2,1))<1||round(line_positions(i,2,1))>sz
                loop = false;
            else
                filaments(round(line_positions(i,1,1)),round(line_positions(i,2,1)),1) = 1;
%                 filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
                if rand()<pr
                    particles.position(j,:) = [line_positions(i,1,1) line_positions(i,2,1)];
                    j = j+1;
                end
            end
        end
    elseif start_points(i,1,2) == 3
        theta = pi*rand();
        while loop
            line_positions(i,1,1) = line_positions(i,1,1) - sin(theta);
            line_positions(i,2,1) = line_positions(i,2,1) + cos(theta);
            if round(line_positions(i,1,1))<1||round(line_positions(i,1,1))>sz||round(line_positions(i,2,1))<1||round(line_positions(i,2,1))>sz
                loop = false;
            else
                filaments(round(line_positions(i,1,1)),round(line_positions(i,2,1)),1) = 1;
%                 filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
                if rand()<pr
                    particles.position(j,:) = [line_positions(i,1,1) line_positions(i,2,1)];
                    j = j+1;
                end
            end
        end
    else
        theta = pi*rand();
        while loop
            line_positions(i,1,1) = line_positions(i,1,1) + cos(theta);
            line_positions(i,2,1) = line_positions(i,2,1) - sin(theta);
            if round(line_positions(i,1,1))<1||round(line_positions(i,1,1))>sz||round(line_positions(i,2,1))<1||round(line_positions(i,2,1))>sz
                loop = false;
            else
                filaments(round(line_positions(i,1,1)),round(line_positions(i,2,1)),1) = 1;
%                 filaments(round(positions(i,1,1)),round(positions(i,2,1)),2) = theta;
                if rand()<pr
                    particles.position(j,:) = [line_positions(i,1,1) line_positions(i,2,1)];
                    j = j+1;
                end
            end
        end
    end
end
particles.position = particles.position(any(particles.position,2),:); % Remove all zero rows

figure()
colormap(pink)
imagesc(filaments)

figure()
image = zeros(sz);
for i = 1:size(particles.position,1)
    image(round(particles.position(i,1)),round(particles.position(i,2))) = 1;
end
imagesc(image)
colormap(pink)
