function [position] = generatePositions(N,sz,w0,varargin)

% default cut-off of non-zero elements in psf
kernelSize = ceil(3*w0);
for ii = 1:2:length(varargin)
    if any(strcmpi(varargin{ii},{'kernelSize'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} > 0
            kernelSize = ceil(varargin{ii+1}*w0);
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
% kernelSize beyond the image edges
position = (1 - kernelSize) + (2*kernelSize+sz-1).*rand(N,2);