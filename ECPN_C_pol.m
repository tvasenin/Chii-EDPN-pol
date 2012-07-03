function P = ECPN_C_pol(rel,E,Wpol)
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
assert(nnz(diag(E))==0, '[ERROR] Assertion failed -- nnz(diag(E))>0 in the start of ECPN_C!');


%if nnz(diag(E)) > 0, disp('HALT'); end
%% 2-3-4 node case

switch n
    case 2 % 2-node graph
        cnt.NUMEL2 = cnt.NUMEL2 + 1;
        %VWpol = Wpol2VWpol_opt(rel,Wpol);
        %P = ECPN_C_numel2_pol(VWpol);
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

%% small precalc
Es = sum(E);
q = sum(Es)/2;
%assert(nnz(rel(Es==2))==0,'out');

%% full/chain/cycle case
%
switch q
    case n*(n-1)/2 % full
       cnt.FULL = cnt.FULL + 1;
       VWpol = Wpol2VWpol_opt(rel,Wpol);
       P = ECPN_full_pol(VWpol);
       return
    case n
        if max(Es)==2 % cycle
            if n>5 % do not highlight small cycles (uncomment for cycles debug!)
                disp(['Found cycle = ', int2str(n)]) 
            end
            cnt.CYCLE = cnt.CYCLE + 1;
            P = ECPN_cycle_pol_v2(rel,E,Wpol);
            return
        end
    case n-1 % tree
        if max(Es)==2 %chain
            disp('[INFO] Chain has been found inside ECPN_C!')
            cnt.CHAIN = cnt.CHAIN + 1;
%            disp(['conncomp = ', int2str(graphconncomp(E,'Directed',false))])
            VWpol = Wpol2VWpol_opt(rel,Wpol);
            P = ECPN_chain_pol(rel,E,VWpol);
            return
        end
end

%

%% (connected) graph with node of n-1 degree
%
if (max(Es) == n-1)  % may be optimized for multiple nodes
    cnt.MAXDEG = cnt.MAXDEG + 1;
    %temporary disabled debug output
    max_mask = Es==(n-1);
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

%% Debug -- checking for hanging trees(or nodes or chains)
%
%also should find trees really good
%clear q
%P = 0;

if ~isempty(find(Es==1,1))
    cnt.HNODES = cnt.HNODES + 1;
%    disp('Found hnodes!');
%    if (q == n-1)
%    cnt.TREE = cnt.TREE + 1;
%        disp(['Found tree = ' int2str(n)])
%    end
%    P = ECPN_hnodes_pol_v3(rel,E,Wpol);
%    return;
    [P, rel, E, Wpol] = ECPN_hnodes_pol_v3_back(rel,E,Wpol);
%    if numel(rel) > 3 %
    if nnz(E) > 0 %still have some edges
        tmp = ECPN_C_pol(rel,E,Wpol);
        P = poly_add(P,tmp);
    end
    return
end

%Es = sum(E); % recalcing
%% Reduction of remaining chains (if any)
%

assert (max(Es) > 2,'[ERROR] Got cycle after cycle check!');
% Assert max(Es) > 2

%detect chains
ind_ch = Es==2;

%have to temporarily disable this optimization
%because of conflicts with reliable nodes reducing
%to enable ASAP

%if nnz(ind_ch) > 0 %have 2 or more nodes with deg = 2
if nnz(E(ind_ch,ind_ch)) % have at least two consecutive nodes of degree 2
    cnt.CHAINRED = cnt.CHAINRED + 1;
    P = ECPN_chainred_pol_v3(rel,E,Wpol);
    return
end

%% Fully reliable connected graph
%
%think about moving to ECPN
assert(nnz(~rel)>0,'[ERROR] Assertion failed -- found fully reliable connected graph in ECPN_C just before branching!');
if nnz(~rel)==0 %V == ones(1,n) %sum(V) == n
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

[~, index] = sort(Es,'descend');

rel = rel(index);
E = E(index,index);
Wpol = Wpol(index,:);

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

%% Trying to find junction points
%
% not sure this will help ATM
% test_ind = true(1,n);
% test_ind(V == 1) = false;
% 
% nodes_test = 1:n;
% nodes_test = wrev(nodes_test(test_ind));
% 
% for i = nodes_test
%     [CompNum ~] = graphconncomp(E([1:i-1 i+1:n],[1:i-1 i+1:n]),'Directed',false);
%     if CompNum > 1
%         red = i;
%         break;
%     end
% end

%%
%
% m is always p (or [1 0]) because we are contracting it!
%m = V(red); % 

assert(~rel(red),'[HALT] Assertion failed: reliable node has been chosen for reduction!');
%% first subgraph
ind = true(1,n);
ind(red) = 0;

%P = P + (1-m) * ECPN(Vred, Ered, Wred);
P = conv2(ECPN_pol(rel(ind),E(ind,ind),Wpol(ind,:)),[-1, 1]);

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
E = ContractNodes(E,red);

%% second subgraph

%P = P + p * ECPN_C(V(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),W(~nodes_rel_ind));
tmp = ECPN_C_pol(rel(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),Wpol(~nodes_rel_ind,:));
P = poly_add(P, [tmp 0]);

%removing nodes_rel;
%V(nodes_rel) = [];
%E(nodes_rel,:) = []; E(:,nodes_rel) = [];
%W(nodes_rel) = [];
%P = P + p * ECPN_C(V,E,W);
%P = (1-m) * ECPN(V(ind), E(ind,ind), W(ind)) + m * ECPN_C(V,E,W);

end