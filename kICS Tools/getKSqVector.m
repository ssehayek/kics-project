function [ksq_sub,iksq_sub] = getKSqVector(r_k,varargin)

for i = 1:length(varargin)
    if strcmpi(varargin{i},{'kSqMin'})
        if isnumeric(varargin{i+1})
            ksq_min = varargin{i+1};
        elseif strcmpi(varargin{i+1},'min')
            % do nothing; min set by default
        else
            warning(['Unknown option for ''',varargin{i},...
                ''', using default options.'])
        end
    elseif strcmpi(varargin{i},{'kSqMax'})
        if isnumeric(varargin{i+1})
            ksq_max = varargin{i+1};
        elseif strcmpi(varargin{i+1},'max')
            % do nothing; max set by default
        else
            warning(['Unknown option for '' ',varargin{i},...
                ''', using default options.'])
        end
    end
end

% compute spatial dimensions of input array
size_y = size(r_k,1);
size_x = size(r_k,2);
%

% get 1-d vector of corresponding lattice points
% treat even & odd cases separately
if mod(size_x,2) == 0
    xgv = 0:size_x/2;
else
    xgv = 0:(size_x-1)/2;
end

if mod(size_y,2) == 0
    ygv = 0:size_y/2;
else
    ygv = 0:(size_y-1)/2;
end
%

% create xy-lattices
[X,Y] = meshgrid(xgv,ygv); 
lattice_sqrd = (X/size_x).^2 + (Y/size_y).^2; % norm squared of lattice
% unique values occuring in lattice sorted into vector
lattice_nums = unique(lattice_sqrd); 

ksq = (2*pi)^2*lattice_nums; % unique |k|^2 values sorted into vector

% lowest index, i, which satisfies kSqVector(i) >= kSqMin
if exist('ksq_min','var')
    iksq_min = find(ksq >= ksq_min,1,'first');
else
    iksq_min = 1;
end
%
% highest index, j, which satisfies kSqVector(j) <= kSqMax
if exist('ksq_max','var')
    iksq_max = find(ksq <= ksq_max,1,'last');
else
    iksq_max = length(ksq);
end
%
ksq_sub = ksq(iksq_min:iksq_max); % all values satisfying kSqMin <= kSqVector <= kSqMax
iksq_sub = iksq_min:iksq_max; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax
