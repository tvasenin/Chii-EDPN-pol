function [ index, c0, cend, c_first, c_last ] = GetIndexedChain( E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

%% assert: no hnodes, not a cycle

n = length(E);
Es = sum(E);
ind_ch = (Es == 2);

E_ch = E;
E_ch(~ind_ch,:) = 0;
E_ch(:,~ind_ch) = 0;
%[CompNum comps] = graphconncomp(E_ch,'Directed',false);
[CompNum comps] = graphalgs('wcc',0,false,E_ch); % shortcut, need to place graphalgs MEX-file to the MATLAB path

c_lens = zeros(0,CompNum);
for i = 1:CompNum
    c_lens(i) = nnz(comps == i);
end
c_lens(c_lens == 1) = 0;

[c_len c_max] = max(c_lens);

%disp([int2str(nnz(c_len)), ' chain(s) found (max length = ', int2str(c_max_len), '): ', int2str(nonzeros(c_len)')]);
c = find(comps == c_max);
if length(c) > 2
    %Es_ch = sum(E_ch);
    %c = c(graphtraverse(E(c,c),find(Es_ch(c)==1,1),'Directed',false)); %traversing chain
    %c = c(graphalgs('dfs',0,false,E(c,c),find(Es_ch(c)==1,1),inf));% shortcut, need to place graphalgs MEX-file to the MATLAB path
    c = c(chaintraverse_fast(E(c,c)));
end
c0   = find(E(c(1),:)   & ~E_ch(c(1),:)  );
cend = find(E(c(end),:) & ~E_ch(c(end),:));
index = 1:n;

%if c0 == cend
%    disp('OLOLO! c0 = cend!');
%end

index([c0 c]) = [];
index = [c0 c index];
c0   = find(index == c0);
cend = find(index == cend);

c_first = 2;
c_last = c_first + c_len - 1;

end

