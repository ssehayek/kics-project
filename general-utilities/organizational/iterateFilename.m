% iterateFilename(filename) appends the lowest unused number to the
% input directory/file string, and returns this newly appended string.
%
% Either absolute, or relative path names may be used, as well as file
% names. Enter directories with trailing slashes.
function [new_filename,n] = iterateFilename(filename,varargin)

% string to prepend before lowest unused number
prep_str = '';
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'prepStr','prependString'}))
        if ischar(varargin{ii+1})
            prep_str = varargin{ii+1};
        else
            warning(['Unknown option for ''',varargin{ii},...
                ''', using default options.'])
        end
    end
end

[pathstr,name,ext] = fileparts(filename);

if isempty(pathstr)
    target_n = @(n) [name,prep_str,num2str(n),ext];
elseif isempty(name) && isempty(ext)
    target_n = @(n) [pathstr,prep_str,num2str(n)];
else
    target_n = @(n) [pathstr,filesep,name,prep_str,num2str(n),ext];
end

n = 1; % n is the lowest number which is not appended to a folder/file
while exist(target_n(n)) ~= 0 % iterate so long as the target directory/file exists
    n = n + 1;
end
new_filename = target_n(n); % target directory with lowest unused n appended
