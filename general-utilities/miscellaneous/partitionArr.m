% splits array, arr, into n_parts, or n_parts-1 equal parts. In the former
% case, if numel(arr) is not divisible by n_parts, then the vector is
% split into n_parts-1 equal parts, and one extra part containing the
% remainder of elements at the end of the vector
%
% INPUT
% arr: vector to be split
% n_parts: request number of parts vector to be partitioned into
% 
% OUTPUT
% arr_parts: cell containing the segmented vector parts
% n_parts: number of parts arr was split into (arr might be split into 
% n_parts-1, instead of n_parts)
%
function [arr_parts,n_parts] = partitionArr(arr,n_parts)

% numel is used to generalize function for an array
switch isscalar(arr)
    case 0 
        % arr is array
        n = numel(arr);
    otherwise
        % if positive scalar integer, n, was input for arr, turn it into an
        % array so that arr(1)=1 and arr(end)=n
        n = arr;
        arr = 1:n;
end
    
r1 = mod(n,n_parts);
switch r1
    % check if n is divisible by n_parts
    case 0 
        % arr is split into n_parts equal parts
        n_rows = n/n_parts;
        % each column of eq_vecs is an equal partition of arr
        eq_vecs = reshape(arr,[n_rows,n_parts]);
        arr_parts = cell(1,n_parts);
        for ii = 1:n_parts
            arr_parts{ii} = eq_vecs(:,ii);
        end
    otherwise
        r2 = mod(n,n_parts-1);
        switch r2
            % check if n is divisible by n_parts-1
            case 0 
                % arr is split into n_parts-1 equal parts
                warning(['array was partitioned into ',num2str(n_parts-1),...
                    ' parts, instead of ',num2str(n_parts),'.'])
                n_parts = n_parts - 1;
                n_rows = n/n_parts;
                eq_vecs = reshape(arr,[n_rows,n_parts]);
                arr_parts = cell(1,n_parts);
                for ii = 1:n_parts
                    arr_parts{ii} = eq_vecs(:,ii);
                end
            otherwise
                % arr is split into n_parts, with last part containing
                % remainder of elements
                %
                % number of elements in each of the n_parts-1 equal
                % partitions (from modulo formula)
                n_rows = (n-r2)/(n_parts-1);
                eq_vecs = reshape(arr(1:n-r2),[n_rows,n_parts-1]);
                arr_parts = cell(1,n_parts);
                for ii = 1:n_parts-1
                    arr_parts{ii} = eq_vecs(:,ii);
                end
                arr_parts{end} = arr(n-r2+1:end);
        end
end