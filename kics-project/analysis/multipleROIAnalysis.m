file = dir('**/ROI*.mat');
filename = {file.name};
filedir = {file.folder};

n = length(filename);
X = zeros(n,4);
for i = 1:n
    file_i = [filedir{i},filesep,filename{i}];
    
    m = matfile(file_i);
    opt_params_i = m.opt_params;
    X(i,:) = opt_params_i;
end