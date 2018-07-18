% overlay Gaussian laser profile onto an image series, I
%
% INPUT
%
% I: image series to overlay laser profile
%
% VARARGIN (if no options are input, this code does nothing)
%
% - choose e^(-2) radius of laser profile
%        'laserWidth' | values: (default Inf) 0 to Inf
% - shift of laser peak in pixels; enter as [X0,Y0]
%        'laserShift' | values: (default [0,0]) any numeric 2-vector
function [J,I,laser_profile] = addLaserProfile(I,varargin)

w = Inf; % default width of Gaussian set
X0 = 0; Y0 = 0;
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'laserWidth','laserRadius'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} >= 0
            w = varargin{ii+1};
        elseif strcmpi(varargin{ii+1},'default')
            % set to default
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    elseif any(strcmpi(varargin{ii},{'laserShift'}))
        if length(varargin{ii+1}) == 2 && all(isnumeric(varargin{ii+1}))
            X0 = varargin{ii+1}(1); Y0 = varargin{ii+1}(2);
        elseif strcmpi(varargin{ii+1},'default')
            % set to default
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    end
end
% center pixel of I
[ctr_pxl_y,ctr_pxl_x] = getCtrPxl(I);
% 
xgv = (1:size(I,2))-ctr_pxl_x;
ygv = (1:size(I,1))-ctr_pxl_y;
% 
[X,Y] = meshgrid(xgv,ygv);

laser_profile = exp(-2*((X-X0).^2+(Y-Y0).^2)/w^2);

J = I.*laser_profile;