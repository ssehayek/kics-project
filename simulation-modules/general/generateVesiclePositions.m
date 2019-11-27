function [dye_positions,vesicles] = generateVesiclePositions(N,sz,w0,varargin)

% default vesicles' radii
radius = 0;
% default number of vesicles
n_vesicles = N;
% default cut-off of non-zero elements in psf
kernelSize = ceil(3*w0);
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'kernelSize'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} > 0
            kernelSize = ceil(varargin{ii+1}*w0);
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    elseif any(strcmpi(varargin{ii},{'vesicleRadius','vesRadius','vesR','rVes'}))
        if isnumeric(varargin{ii+1})
            radius = varargin{ii+1};
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    elseif any(strcmpi(varargin{ii},{'nVesicles','numVesicles','numVes','nVes'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} > 0
            n_vesicles = varargin{ii+1};
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    elseif any(strcmpi(varargin{ii},{'boundCond','kerType','kernelType'}))
        % do nothing; avoids warning for later use
    else
        warning(['unknown varargin input ''',varargin{ii},'''.'])
    end
end
% generate random positions of dyes, allowing for positions to be at most
% ( kernelSize - radius ) beyond the image edges
vesicles.positions = (1 - kernelSize + radius) + ...
    (2*(kernelSize-radius)+sz-1).*rand(n_vesicles,2);
% randomly associate dyes to vesicles
vesicles.nDyes = mnrnd(N,ones(1,n_vesicles)*1/n_vesicles);
% cumulative dye count over indexed vesicles
dye_csum = [0,cumsum(vesicles.nDyes)];
% generate dye positions relative to vesicle centers
dye_radius = radius*sqrt(rand([N,1]));
dye_theta = 2*pi*rand([N,1]);

% absolute dye positions
vesicles.dyePositions = cell(1,n_vesicles);
dye_positions = zeros(N,2);
for m = 1:n_vesicles
    m1 = dye_csum(m) + 1;
    m2 = dye_csum(m+1);
    % 
    vesicles.dyePositions{m} = vesicles.positions(m,:) + ...
        [dye_radius(m1:m2).*cos(dye_theta(m1:m2)),...
        dye_radius(m1:m2).*sin(dye_theta(m1:m2))];
    dye_positions(m1:m2,:) = vesicles.dyePositions{m};
end