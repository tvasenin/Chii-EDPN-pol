function [ lmat, pmat, mcnt ] = GetHnodesQueue( E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% naive approach, can be rewritten to be much faster

n = length(E);
Es = sum(E);
q = sum(Es)/2;
lmat = zeros(n-1,n-1); % n-1 nodes max!
pmat = zeros(n-1,n-1);
mcnt = zeros(1,n-1);
iter = 0;

ileaves = (Es==1);

while nnz(ileaves)
    iter = iter + 1;

    pars = find(sum(E(ileaves,:),1)); % explicit sum by dim 1 !
    
    if (q == 1) % only 2 nodes left, both are leaves
        % we have a tree, congratulations
        pars = pars(1);
    end
    
    m = length(pars);
    
    leaves2 = zeros(1,m);
    for i = 1:m
        leaves2(i) = find(E(pars(i),:) & ileaves,1);
    end

    mcnt(iter)      = m;
    lmat(iter, 1:m) = leaves2;
    pmat(iter, 1:m) = pars;
    
    Es(leaves2) = 0;
    Es(pars) = Es(pars) - 1;

    q = q - m;

    E(leaves2,:) = 0;
    E(:,leaves2) = 0;

    ileaves = (Es==1);
end

% truncate output
mcnt = mcnt(1:iter);
lmat = lmat(1:iter, 1:max(mcnt));
pmat = pmat(1:iter, 1:max(mcnt));

end
