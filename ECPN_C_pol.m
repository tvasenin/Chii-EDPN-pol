function P = ECPN_C_pol(rel,E,Wpol, cut_idx)
%ECPN_C Calculates ECPN for given connected(!) G(V,E,W)
%   Assume that G is connected and has at least 1 edge!

global cnt

%global HIT MISS

cnt.TOTAL = cnt.TOTAL + 1;
if mod(cnt.TOTAL,5000) == 0
    disp(['TOTAL  cycles passed: ',int2str(cnt.TOTAL)]);
    disp(['BRANCH cycles passed: ',int2str(cnt.BRANCHING)]);
end

%% small precalc
n = length(rel);
assert(~any(diag(E)), '[ERROR] Assertion failed -- nnz(diag(E))>0 in the start of ECPN_C!');

%% 2-3-4 node case

switch n
    case 2 % 2-node graph
        cnt.NUMEL2 = cnt.NUMEL2 + 1;
        P = ECPN_C_numel2_pol_v2(rel,Wpol);
        return
    case 3 % 3-node graph
        cnt.NUMEL3 = cnt.NUMEL3 + 1;
        VWpol = Wpol2VWpol_opt(rel,Wpol);
        P = ECPN_C_numel3_pol(rel,E,VWpol);
        return
    case 4 % 4-node graph
        cnt.NUMEL4 = cnt.NUMEL4 + 1;
        P = ECPN_C_numel4_pol(rel,E,Wpol);
        return
end

%% now we can assume that we have connected graph with 5 or more nodes

%% precalc
Es = sum(E);
%assert(~any(rel(Es==2)),'out');

%% (connected) graph with node of n-1 degree, but not full
%
if (max(Es) == n-1) && any(Es~=(n-1)) % should not be full, may be optimized for multiple nodes
    cnt.MAXDEG = cnt.MAXDEG + 1;
    %temporary disabled debug output
    %max_mask = Es==(n-1);
    %disp(['[INFO] found ',int2str(nnz(max_mask)),' nodes of degree n-1!'])

    VWpol = Wpol2VWpol(rel,Wpol);
    
    if true %nnz(max_mask) == 1 % single node
        maxdeg = find(Es==(n-1),1);
        mask = [1:maxdeg-1 maxdeg+1:n];
     
        %m = V(maxdeg)
        %Wpol  = sym2poly_arr(W);
        %VWpol = sym2poly_arr(VW); can't do this because of unpredictable result length!
        %VWpol = Wpol2VWpol(rel,Wpol);
        %VW(maxdeg)  = W(maxdeg); % VW will be used only if maxdeg-th node works;
        VWpol(maxdeg,:) = [0 Wpol(maxdeg,:)];
        
        P = ECPN_full_pol(VWpol);
        if  ~rel(maxdeg) %m == p
            P = [P 0]; %P = p * ECPN_full(VW);
            tmp = ECPN_pol(rel(mask),E(mask,mask),Wpol(mask,:));
            %tmp = conv2([-1,1],tmp);
            tmp = -[tmp 0] + [0 tmp];
            P = poly_add(P,tmp);
        end
            
        %P = m * feval(symengine, 'ECPN_full_mu',VW) + (1-m) * ECPN(V(mask),E(mask,mask),W(mask));
    else %multiple nodes
        %Qp = prod(1-V(max_mask));
        %Ec = ContractNodes(E,find(max_mask));   
        %P =        Qp  * ECPN(V(~max_mask), E(~max_mask,~max_mask),W(~max_mask));
        %P = P + (1-Qp) * ECPN(V(~max_mask),Ec(~max_mask,~max_mask),W(~max_mask));
        %P = P + feval(symengine, 'ECPN_full_mu',VW(max_mask));
        %P = P + sum(VW(max_mask)) * sum(VW(~max_mask));
    end
    return
end

%% precalc
q = sum(Es)/2;

%% full/cycle/chain/tree case
%
switch q
    case n*(n-1)/2 % full
       cnt.FULL = cnt.FULL + 1;
       VWpol = Wpol2VWpol_opt(rel,Wpol);
       P = ECPN_full_pol(VWpol);
       return
    case n %cycle with possible hnodes
        cnt.CYCLE = cnt.CYCLE + 1;
        if max(Es)>2 % cycle with hnodes
            cnt.HNODES = cnt.HNODES + 1;
            [P, rel, E, Wpol] = ECPN_hnodes_pol_v4_back(rel,E,Wpol);
            n = numel(rel);
            P = poly_add(P,ECPN_cycle_pol_v2(rel,E,Wpol));
        else % max(Es)==2 % cycle
            P = ECPN_cycle_pol_v2(rel,E,Wpol);
        end
        if n>5 % do not highlight small cycles (uncomment for cycles debug!)
            disp(['Found cycle = ', int2str(n)]) 
        end
        return
    case n-1 % tree
        if max(Es)==2 %chain
            disp('[INFO] Chain has been found inside ECPN_C!')
            cnt.CHAIN = cnt.CHAIN + 1;
%            disp(['conncomp = ', int2str(graphconncomp(E,'Directed',false))])
            VWpol = Wpol2VWpol_opt(rel,Wpol);
            P = ECPN_chain_pol(rel,E,VWpol);
        else %tree and not chain
            disp('[INFO] Tree has been found inside ECPN_C!')
            cnt.TREE = cnt.TREE + 1;
            [P, ~, ~, ~] = ECPN_hnodes_pol_v4_back(rel,E,Wpol);
        end
        return        
end

%% Debug -- checking for hanging trees(or nodes or chains)
%
%also should find trees really good
%clear q
%P = 0;

if ~isempty(find(Es==1,1))
    cnt.HNODES = cnt.HNODES + 1;
    [P, rel, E, Wpol] = ECPN_hnodes_pol_v4_back(rel,E,Wpol);
%    if numel(rel) > 3 %
    if find(E,1) %still have some edges
        tmp = ECPN_C_pol(rel,E,Wpol);
        P = poly_add(P,tmp);
    end
    return
end

%Es = sum(E); % recalcing
%% Reduction of remaining chains (if any)
%

assert (max(Es) > 2,'[ERROR] Got cycle after cycle check!');

%detect chains
ind_ch = Es==2;

%have to temporarily disable this optimization
%because of conflicts with reliable nodes reducing
%to enable ASAP

%if any(ind_ch) %have 2 or more nodes with deg = 2
if find(E(ind_ch,ind_ch),1) % have at least two consecutive nodes of degree 2
    cnt.CHAINRED = cnt.CHAINRED + 1;
    P = ECPN_chainred_pol_v3(rel,E,Wpol);
    return
end

%% Fully reliable connected graph
%
%think about moving to ECPN
assert(~all(rel),'[ERROR] Assertion failed -- found fully reliable connected graph in ECPN_C just before branching!');
if all(rel) %V == ones(1,n) %sum(V) == n
    cnt.RELIABLE = cnt.RELIABLE + 1;
    P = ECPN_full_pol(Wpol); %using W instead of WV
    return
end

%%
%
cnt.BRANCHING = cnt.BRANCHING + 1;


%% Sorting nodes by degree;
%

%[Es, index] = sort(Es,'descend');
%index = colperm(E);
%Es = Es(index);

%[Es, index] = sort(Es,'descend');
[~, index] = sort(Es,'descend');

rel = rel(index);
E = E(index,index);
Wpol = Wpol(index,:);

if nargin < 4
    cut_idx = VertexCuts(E);
else
    cut_idx = cut_idx(index);
end

%ascending order slows up things!
%needs testing

%% choosing node for reduction
%

red = find(~rel,1); % but maybe it's better to find another high-degree node

% %if false
% %if (q - n) == 0 % need to remove   1 edge  to get a tree
%                  % never het here because no hanging nodes
% if (q - n) == 1 % need to remove 1-2 edges to get a tree
% %if (q - n) <= 2 % need to remove 1-3 edges to get a tree
% %if (q - n) <= 3 % need to remove 1-4 edges to get a tree
%     [sptree] = graphminspantree(E,1);
%     Esdiff = Es - ( sum(sptree) + sum(sptree,2)' );
%     red1 = find(Esdiff & (V~=1),1);
%     if isempty(red1) || (red1 == red)
%         MISS = MISS + 1;
%         disp(['fap-ololo ', int2str(q-n), ' MISSED']);
%     else
%         HIT  =  HIT + 1;
%         disp(['fap-ololo ', int2str(q-n)]);
%         red = red1;
%     end
% end

%% if we have one node with degree n-2, we can reduce by another
%tmp = find(Es == n-2,1); %will be used until splitting to ECPN and ECPN_C

% if false %~isempty(tmp) %wait for implementation of contractions;
%if false
% if max(Es) == n-2
%     tmp2 = E(1,:); %tmp = 1 because Es is sorted
%     tmp2(1) = 1;
%     tmp3 = find((tmp2 == 0) & (V~=1));
%     if ~isempty(tmp3)
%         red = tmp3;
%         disp('fap_good');
%     else disp('fap_bad');
%     end
% end

%% Trying to find best junction (unreliable) node
%

junc_cand = find((~rel) & cut_idx);

if ~isempty(junc_cand)
    comp_weight(1:n) = (n-1)^2;
    for i = junc_cand
        [~, c_lens] = components_mex(E([1:i-1 i+1:n],[1:i-1 i+1:n]));
        if (length(c_lens) > 1)
            comp_weight(i) = sum(c_lens.^2);
        end
    end
    %comp_weight = comp_weight / (n-1)^2; % assert: n > 1 :)
    [~,red] = min(comp_weight);
%    disp('[INFO] Found junction point!');
end

%%
%
% m is always p (or [1 0]) because we are contracting it!
%m = V(red); % 

assert(~rel(red),'[HALT] Assertion failed: reliable node has been chosen for reduction!');
%% first subgraph (failed node)
ind = true(1,n);
ind(red) = 0;

%P = P + (1-m) * ECPN(Vred, Ered, Wred);
if ~cut_idx(red)% red is not cut, so after removing red graph will stay connected
    P = ECPN_C_pol(rel(ind),E(ind,ind),Wpol(ind,:)); % may get new cut nodes, can't pass cut_idx! 
else            % red is cut, so we definitely have multigraph
    if any(cut_idx & E(red,:)) %there are some cut nodes adjacent to red
        P = ECPN_pol(rel(ind),E(ind,ind),Wpol(ind,:)); % may lose some cut nodes (red neighbours)!
    else                       %there are some cut nodes adjacent to red
        P = ECPN_pol(rel(ind),E(ind,ind),Wpol(ind,:), cut_idx(ind)); % red is cut node and there are no adjacent cut nodes
    end
end
P = conv2(P,[-1, 1]);

%% searching for another adjacent reliable nodes
%  
%V(red) = 1;
rel(red) = 1;
nodes_rel_ind = rel & E(red,:);
nodes_rel = find(nodes_rel_ind);


if ~isempty(nodes_rel)
    E = ContractNodes(E,nodes_rel);
    Wtmp_pol = Wpol([nodes_rel red],:);
    
    tmp = [ECPN_full_pol(Wtmp_pol) 0];
    P = poly_add(P,tmp);
    
    %W(red) = W(red) + sum(W(nodes_rel));
    Wpol(red,:) = Wpol(red,:) + sum(Wpol(nodes_rel,:),1); % !!! CAREFUL !!!
end

%% contracting E without deleting red node
E = ContractNodes(E,red); % maybe we can contract red with nodes_rel

%% second subgraph (reliable node)

%Will not get new cuts after contracting,
%but need strict check which nodes will ctop being cuts after contracting
%before using cut_idx here!
%P = P + p * ECPN_C(V(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),W(~nodes_rel_ind));
tmp = ECPN_C_pol(rel(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),Wpol(~nodes_rel_ind,:));
P = poly_add(P, [tmp 0]);

end