% created by Simon Sehayek
% credit to Hugo B. Brandao
%
% this function outputs the psf kernel of a molecule located at pos
% (format: [pos_y,pos_x]), in an image series with spatial dimensions
% specified by input variable J, with e-2 radius of w0 (assumes gaussian psf).
%
% INPUT VARIABLES
%
%   J: image series
%   pos: position of molecule (format: [pos_y,pos_x])
%   w0: Gaussian PSF e-2 radius
%
% OUTPUT VARIABLES
%
%   img_kernel: PSF kernel of molecule with position pos
%
% VARARGIN (NAME/VALUE PAIRS)
%
% 'boundCond'
% -------------------------------------------------------------------------
% 'hard' (default) | 'periodic'
% -------------------------------------------------------------------------
% Specifies boundary condition at edges. 'hard' means the PSF does not go
% beyond the edges of the image; 'periodic' means the part of the PSF
% exceeding the image boundaries appears on the opposite side of the image
%
% 'kernelSize'
% -------------------------------------------------------------------------
% ceil(3*w0) (default) | cell-array
% -------------------------------------------------------------------------
% Specifies cut-off for PSF in kernel (in units of w0)
%
% 'kerType'
% -------------------------------------------------------------------------
% 'integrate' (default) | 'legacy'
% -------------------------------------------------------------------------
% Specifies PSF. 'integrate' means the relative intensity value in each
% pixel is the result of integrating the PSF along pixel widths; 'legacy'
% means the relative intensity value in each pixel is taken to be the value
% of the PSF at the center of the respective pixel
%
function [img_kernel] = getImgKernel(J,pos,w0,varargin)

% default cut-off of non-zero elements in psf
kernelSize = ceil(3*w0);
% default boundary condition
bound_cond = 'hard';
% choose to either evaluate psf at a point ('legacy' option) or integrate
% over psf over pixels ('integrate')
ker_type = 'integrate';
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'kernelSize'}))
        if isnumeric(varargin{ii+1}) && varargin{ii+1} > 0
            kernelSize = ceil(varargin{ii+1}*w0);
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    elseif any(strcmpi(varargin{ii},{'boundCond'}))
        if any(strcmpi(varargin{ii+1},{'periodic'}))
            bound_cond = 'periodic';
        elseif any(strcmpi(varargin{ii+1},{'hard'}))
            % default
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    elseif any(strcmpi(varargin{ii},{'kerType','kernelType'}))
        if any(strcmpi(varargin{ii+1},{'legacy'}))
            ker_type = 'legacy';
        elseif any(strcmpi(varargin{ii+1},{'integrate'}))
            % default
        else
            warning(['invalid option for varargin: ',varargin{ii}]);
        end
    end
end

% relevant (i.e. non-zero) coordinates relative to psf center
[x,y] = meshgrid(-kernelSize:kernelSize);

% spatial image series sizes
size_y = size(J,1); 
size_x = size(J,2);

% molecule's position
pos_y = pos(1);
pos_x = pos(2);

% position relative to center
dy = -round(pos_y)+pos_y; 
dx = -round(pos_x)+pos_x; 

if strcmp(ker_type,'legacy')
    % gaussian psf form (intensity amplitude set to 1)
    arg = -2*((x-dx).^2 + (y-dy).^2)/w0^2;
    psf_kernel = exp(arg);
else
    % integrate over gaussian psf (normalized to area 1)
    %
    % x-integral
    x_int = 1/sqrt(8)*(erf(sqrt(2)*(x-1/2-dx)/w0)-erf(sqrt(2)*(x+1/2-dx)/w0));
    % y-integral
    y_int = 1/sqrt(8)*(erf(sqrt(2)*(y-1/2-dy)/w0)-erf(sqrt(2)*(y+1/2-dy)/w0));
    % integrated psf; 2 is normalization factor
    psf_kernel = 2*x_int.*y_int;
end

% rounded coordinates of PSF
ycoor = y + round(pos_y);
xcoor = x + round(pos_x);

% set negligible elements of PSF profile to 0
nonZeroEl = find(psf_kernel);

psf_kernel = psf_kernel(nonZeroEl);
ycoor = ycoor(nonZeroEl);
xcoor = xcoor(nonZeroEl);

if strcmp(bound_cond,'periodic')
    % get coordinates subject to periodic boundary conditions
    %
    %     ycoor = mod(y(nonZeroEl) + round(pos_y),size_y);
    ycoor = mod(ycoor,size_y);
    %     xcoor = mod(x(nonZeroEl) + round(pos_x),size_x);
    xcoor = mod(xcoor,size_x);
    % for MATLAB indexing convention
    xcoor(xcoor==0) = size_x;
    ycoor(ycoor==0) = size_y;
else
    % find kernel positions which exceed boundary and crop them out
    %
    ycoor_out = (ycoor < 1) | (ycoor > size_y);
    xcoor_out = (xcoor < 1) | (xcoor > size_x);
    % coordinates to crop out
    crop_coors = xcoor_out | ycoor_out;
    % cropping
    ycoor(crop_coors) = [];
    xcoor(crop_coors) = [];
    psf_kernel(crop_coors) = [];
end

% embed kernel in rest of image
% img_kernel = full(sparse(ycoor,xcoor,psf_kernel(nonZeroEl),size_y,size_x));
img_kernel = full(sparse(ycoor,xcoor,psf_kernel,size_y,size_x));
