% iterateFilename(filename) appends the lowest unused number to the 
% input directory/file string, and returns this newly appended string.
%
% either absolute, or relative path names may be used, as well as file
% names.
function [new_filename,n] = iterateFilename(filename)

[pathstr,name,ext] = fileparts(filename);
    
if isempty(pathstr)
    target_n = @(n) [name,num2str(n),ext];
else
    target_n = @(n) [pathstr,filesep,name,num2str(n),ext];
end

n = 1; % n is the lowest number which is not appended to a folder/file
while exist(target_n(n)) ~= 0 % iterate so long as the target directory/file exists
    n = n + 1;
end
new_filename = target_n(n); % target directory with lowest n appended 
