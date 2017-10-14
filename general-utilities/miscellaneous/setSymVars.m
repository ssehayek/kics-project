% setSymVars(...) returns structure s with field names given by var_names,
% with symbolic variables under field names given by sym_vars, and nunmeric
% values under remaining field names. Also returns symbolic parameter names
% and corresponding values in cell arrays sym_vars and sym_vals.
%
% INPUT PARAMS
% var_names: cell containing all variable names
% all_vals: cell containing parameter values, ordered according to
%   var_names
% sym_vars: cell subset of var_names which will have symbolic variables
%   under their field names

%-1 function [s,sym_vars,sym_vals] = setSymVars(var_names,all_vals,sym_vars)
function [s] = setSymVars(var_names,all_vals,sym_vars,varargin)

use_vpa = 0;
for ii = 1:length(varargin)
    if any(strcmpi(varargin{ii},{'useVPA'}))
        if varargin{ii+1} == 1 
            use_vpa = 1;
        end
    end
end

% intersection of var_names and sym_vars, in case sym_vars is not a subset
% of var_names
[sym_vars,i_sym] = intersect(var_names,sym_vars);
sym_vals = {all_vals{i_sym}};

% remaining non-symbolic variables
[nonsym_vars,i_nonsym] = setdiff(var_names,sym_vars);
nonsym_vals = {all_vals{i_nonsym}};

% variables not in var_names throw warning
non_vars = setdiff(sym_vars,var_names);
non_vars = strjoin(non_vars,', ');
if ~isempty(non_vars)
    warning(['unknown variables: ',non_vars,' passed in sym_vars'])
end

% store symbolic variables in s
for ii = 1:length(sym_vars)
    switch use_vpa
        case 0
            s.(sym_vars{ii}) = sym(sym_vals{ii});
        otherwise
            s.(sym_vars{ii}) = vpa(sym_vals{ii});
    end
end
% store non-symbolic variables in s
for ii = 1:length(nonsym_vars)
    s.(nonsym_vars{ii}) = nonsym_vals{ii};
end