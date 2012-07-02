function [b,ndx] = unique_fast(a)
%UNIQUE_FAST Set unique [fast for NeXT].
%   B = UNIQUE_FAST(A) for the array A returns the same values as in A
%   but with no repetitions. B will also be sorted.

numelA = numel(a);
        
% Handle empty: no elements.
if (numelA == 0)
    % Predefine b to be of the correct type.
    b = a([]);
    if max(size(a)) > 0
        b = reshape(b,0,1);
        ndx = zeros(0,1);
    else
        ndx = [];
    end
elseif (numelA == 1)
    % Scalar A: return the existing value of A.
    b = a;
    ndx = 1;
    % General handling.
else
    % Convert to columns
    a = a(:);
    % Sort
    [b,ndx] = sort(a);
    % d indicates the location of non-matching entries.
    d = diff(b) ~= 0;
    d(numelA,1) = true; % Final element is always a member of unique list.
    b = b(d);         % Create unique list by indexing into sorted list.
    % Create indices
    ndx = ndx(d);
end
end
