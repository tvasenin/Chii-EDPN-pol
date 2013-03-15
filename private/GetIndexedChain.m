function [ index, c0, cend, c_first, c_last ] = GetIndexedChain( E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% assert: no hnodes, not a cycle

n = length(E);
ind_ch = (sum(E) == 2);

E_ch = E;
E_ch(~ind_ch,:) = 0;
E_ch(:,~ind_ch) = 0;
[comps c_lens] = components_mex(E_ch); % shortcut, need to place MEX-file to the MATLAB path

[c_len c_max] = max(c_lens);

assert(c_len > 1, 'c_len should be greater than 1!');

%disp([int2str(nnz(c_len)), ' chain(s) found (max length = ', int2str(c_max_len), '): ', int2str(nonzeros(c_len)')]);
c_ind = (comps' == c_max);
c = find(c_ind,c_len);
if (c_len > 2) % need to traverse chain!
    Ecc = E(c_ind,c_ind);
    [~, index2] = sort(dfs_mex(Ecc,find(sum(Ecc)==1,1),0,0)); % shortcut, need to place MEX-file to the MATLAB path
    c = c(index2);
end

c0   = find(E(c(1)  ,:) & ~c_ind, 1);
cend = find(E(c(end),:) & ~c_ind, 1);

%if c0 == cend
%    disp('OLOLO! c0 = cend!');
%end

index = 1:n;
index([c0 c]) = [];

index = [c0 c index];
c0   = find(index == c0);
cend = find(index == cend);

c_first = 2;
c_last = c_first + c_len - 1;

end

