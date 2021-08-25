% Written by Jason S. Leith and Alexander Verge. Additions by Simon
% Sehayek.
function [ stack, iFrame ] = im3read(filename,varargin)
%   Reads frames (vector) images from filename via imread.
%  im3read historically (pre-2015) converted to double.  Now allows
%  preserving orig data type with param 'numtype' set to 'orig'.
%
% im3read(filename) reads up to the first 1000 frames of filename as a double.
%
% im3read(filename,A), where A is 1D numeric array, reads frames A of
% filename.
%
% im3read(filename,A,'ROI',roi) returns the movie in a specific region
% specified by roi. roi is a 2x2 array with the first row containing the
% minimum and maximum indices of interest, respectively, in the horizontal direction and
% the second row being for the vertical direction.
%
% im3read(filename,Param,Value) reads up to the first 1000 frames of
% filename and applies Param,Value inputs:
%   Verbosity
%   TiffWarnings
%   NumType - can be a numeric class or 'orig'.  Reads the image as the
%   specified numeric type.  If set to 'orig' uses the image file's numeric
%   type.  Default is 'double' -- i.e., im3read will convert to double
%   *unless* a param,value pair 'NumType','orig' (or 'uint8', or 'uint16',
%   etc.) is used.

defaultFTR = 1:1e4;
if nargin < 2
    framesToRead = defaultFTR; %increase to increase max frame w/ no specification of length
elseif isnumeric(varargin{1})
    framesToRead = varargin{1};
else
    framesToRead = defaultFTR;
end

verbosity = 1;
tiffWarnings = 0;
numType = 'double'; % this was the implicit default before 2015, %
% since fcn 'zeros' creates matrix of double zeros
% by default, so keeping it that way.
roiBool = 0;
for iVA = 1:length(varargin)
    if strcmpi(varargin{iVA}, 'Verbosity')
        if strcmpi(varargin{iVA+1}, 'quiet')
            verbosity = 0;
        else
            verbosity = 1;
        end
    elseif any(strcmpi(varargin(iVA), {'TiffWarning', 'TiffWarnings'}))
        tiffWarnings = parseBoolean(varargin{iVA+1});
    elseif any(strcmpi(varargin(iVA), {'NumType', 'NumericType', 'NumericClass', 'DataType'}))
        numType = varargin{iVA+1};
    elseif any(strcmpi(varargin{iVA},{'regionOfInterest','ROI','selectSize','size','region'}))
        roi = varargin{iVA+1};
        roi_x = roi(1,:);
        roi_y = roi(2,:);
        roiBool = 1;
    end
end
if tiffWarnings
    % do nothing
else
    warning('off', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
end
varargin = {};

firstFrame = imread(filename);
dataInfo = whos('firstFrame');
if any(strcmpi(numType, {'orig','original'}))
    numType = dataInfo.class;
else
    %keep either default (double) or specified numeric type
end
if roiBool == 0
    stack = zeros([dataInfo.size length(framesToRead)],numType);
else
    stack = zeros(roi_y(2)-roi_y(1)+1,roi_x(2)-roi_x(1)+1,length(framesToRead),numType);
end

for iFrame = 1:length(framesToRead)
    try
        if roiBool == 0
            stack(:,:,iFrame) = imread(filename, 'Index', framesToRead(iFrame), varargin{:}); %if varargin is made empty, why pass it on?
        else
            tempIm = imread(filename, 'Index', framesToRead(iFrame), varargin{:});
            stack(:,:,iFrame) = tempIm(roi_y(1):roi_y(2),roi_x(1):roi_x(2));
            clear tempIm
        end
    catch err
        if verbosity
            disp(['error caught: ' err.identifier]);
        end
        if strcmp(err.identifier, 'MATLAB:imagesci:rtifc:invalidDirIndex')
            if verbosity
                disp(['Terminating inreading before invalid index ' int2str(framesToRead(iFrame)) '.']);
            end
            iFrame = iFrame - 1; % This is fine because you'll break out of the FOR loop before getting
            % to the top of the loop where iFrame would be reset to the next outer-loop value.
            break;
            
            %stack(:,:,iFrame) = zeros([y, x, 1]);
            %disp(['Frame with invalid index ' int2str(framesToRead(iFrame)) ' filled with zeros.']);
        elseif strcmp(err.identifier, 'MATLAB:imagesci:tiffmexutils:libtiffError')
            if verbosity
                disp(['Terminating inreading before final index ' int2str(framesToRead(iFrame)) ' due to "zero tag directory" malarkey.']);
            end
            iFrame = iFrame - 1;
            break;
        else
            disp('error rethrown')
            rethrow(err)
        end
    end
end

%iFrame will be last successfully added frame.
stack = stack(:,:,1:iFrame);
if ~tiffWarnings
    warning('on', 'MATLAB:imagesci:tiffmexutils:libtiffErrorAsWarning');
end

end
