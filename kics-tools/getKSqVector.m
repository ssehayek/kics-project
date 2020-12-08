function [ksq_sub,iksq_sub,lattice_inds] = getKSqVector(r_k,varargin)

ksq_min = 'min';
ksq_max = 'max';
use_zero = 0;
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
    elseif any(strcmpi(varargin{i},{'useZero','includeZero'}))
        if isnumeric(varargin{i+1}) && any(varargin{i+1} == [0,1])
            use_zero = varargin{i+1};
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
    xgv = -size_x/2:size_x/2-1;
else
    xgv = -(size_x-1)/2:(size_x-1)/2;
end

if mod(size_x,2) == 0
    ygv = -size_y/2:size_y/2-1;
else
    ygv = -(size_y-1)/2:(size_y-1)/2;
end
%

% create xy-lattices
[X,Y] = meshgrid(xgv,ygv);
% norm squared of lattice
lattice_sqrd = (2*pi)^2*((X/size_x).^2 + (Y/size_y).^2);
% unique values occuring in lattice sorted into vector
ksq = unique(lattice_sqrd);

% lowest index, i, which satisfies kSqVector(i) >= kSqMin
if strcmpi(ksq_min,'min')
    iksq_min = 1;
    ksq_min = ksq(1);
else
    iksq_min = find(ksq >= ksq_min,1,'first');
end
%
% highest index, j, which satisfies kSqVector(j) <= kSqMax
if strcmpi(ksq_max,'max')
    iksq_max = length(ksq);
    ksq_max = ksq(end);
else
    iksq_max = find(ksq <= ksq_max,1,'last');
end
%
ksq_sub = ksq(iksq_min:iksq_max); % all values satisfying kSqMin <= kSqVector <= kSqMax
iksq_sub = iksq_min:iksq_max; % all indices satisfying kSqMin <= kSqVector(i) <= kSqMax

if nargout == 3
    % get corresponding lattice indices in the range [ksq_min,ksq_max]
    lattice_inds = find( lattice_sqrd >= ksq_min & lattice_sqrd <= ksq_max );
end

if use_zero == 0 && ksq_sub(1) == 0
    ksq_sub(1) = [];
    iksq_sub(1) = [];
    if nargout == 3
        lattice_inds(1) = [];
    end
end