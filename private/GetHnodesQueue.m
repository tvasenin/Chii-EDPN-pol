function [ lmat, pmat, mcnt ] = GetHnodesQueue( E )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

% naive approach, can be rewritten to be much faster

n = length(E);
Es = full(sum(E));
iter = 0;

ileaves = (Es==1);
nnzi = nnz(ileaves);

while nnzi
    iter = iter + 1;
  
    switch nnzi
        case 1
            leaves2 = find(ileaves,1);
            pars = find(E(leaves2,:),1);
            m = 1;
        case 2
            if (nnz(E) == 2) % only 2 nodes left, both are leaves
            % we have a tree, congratulations
                leaves2 = find(ileaves,1);
                pars = find(E(leaves2,:),1);
                m = 1;
            else
                pars = find(any(E(ileaves,:),1)); % explicitly by dim 1 !
                m = length(pars);
%               if (nnz(ileaves) ~= m) % need to search back for relevant leaves
                leaves2 = zeros(1,m);
                for i = 1:m
                    leaves2(i) = find(E(pars(i),:) & ileaves,1);
                end
            end
        otherwise
            pars = find(any(E(ileaves,:),1)); % explicitly by dim 1 !
            m = length(pars);
%           if (nnz(ileaves) ~= m) % need to search back for relevant leaves
            leaves2 = zeros(1,m);
            for i = 1:m
                leaves2(i) = find(E(pars(i),:) & ileaves,1);
            end
    end

    if (iter == 1)
        %initialize lmat and pmat
        mcnt = zeros(1,n-m); % n-m iters max
        lmat = zeros(m,n-m);
        pmat = zeros(m,n-m);
    end

    mcnt(iter)      = m;
    lmat(iter, 1:m) = leaves2;
    pmat(iter, 1:m) = pars;

    Es(leaves2) = 0;
    Es(pars) = Es(pars) - 1;

    E(leaves2,:) = 0;
    E(:,leaves2) = 0;

    ileaves = (Es==1);
    nnzi = nnz(ileaves);
end

% truncate output
mcnt = mcnt(1:iter);
lmat = lmat(1:iter,:);
pmat = pmat(1:iter,:);

end
