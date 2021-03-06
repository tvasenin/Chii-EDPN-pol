function P = ECPN_C(V,E,W)
%ECPN_C Calculates ECPN for given connected(!) G(V,E,W)
%   Assume that G is connected and has >1 nodes!

global cntNUMEL2 cntNUMEL3 cntNUMEL4 cntFULL cntCHAIN cntCYCLE ...
    cntMAXDEG cntTREE cntHNODES cntRELIABLE cntBRANCHING cntTOTAL

%global HIT MISS

cntTOTAL = cntTOTAL + 1;
if mod(cntTOTAL,500) == 0
    disp(['TOTAL  cycles passed: ',int2str(cntTOTAL)]);
    disp(['BRANCH cycles passed: ',int2str(cntBRANCHING)]);
    
end

%% small precalc
syms p
rel = (V == 1);
n = length(E);
Wpol = sym2poly_arr(W);
VWpol = Wpol2VWpol(rel,Wpol);
assert(nnz(diag(E))==0, '[ERROR] Assertion failed -- nnz(diag(E))>0 in the start of ECPN_C!');
%if nnz(diag(E)) > 0, disp('HALT'); end
%% 2-3-4 node case

switch n
    case 2 % 2-node graph
        cntNUMEL2 = cntNUMEL2 + 1;
        P = ECPN_C_numel2_pol(VWpol);
        P = poly2sym (P,p);
        return
    case 3 % 3-node graph
        cntNUMEL3 = cntNUMEL3 + 1;
        P = ECPN_C_numel3_pol(rel,E,VWpol);
        P = poly2sym (P,p);
        return
    case 4 % 4-node graph
        cntNUMEL4 = cntNUMEL4 + 1;
        P = ECPN_C_numel4_pol(rel,E,Wpol);
        P = poly2sym (P,p);
        return
end

%% now we can assume that we have connected graph with 5 or more nodes

%% small precalc
Es = sum(E);
q = sum(Es)/2;
        

%% full/chain/cycle case
%
switch q
    case n*(n-1)/2 % full
       cntFULL = cntFULL + 1;
       P = ECPN_full_pol(VWpol);
       P = poly2sym (P,p);
       return
    case n
        if max(Es)==2 % cycle
            if n>5 % do not highlight small cycles (uncomment for cycles debug!)
                disp(['Found cycle = ', int2str(n)]) 
            end
            cntCYCLE = cntCYCLE + 1;
            P = ECPN_cycle_pol(rel,E,Wpol);
            P = poly2sym (P,p);
            return
        end
    case n-1 % tree
        if max(Es)==2 %chain
            disp('[INFO] Chain has been found inside ECPN_C!')
            cntCHAIN = cntCHAIN + 1;
%            disp(['conncomp = ', int2str(graphconncomp(E,'Directed',false))])
            P = ECPN_chain_pol(rel,E,VWpol);
            P = poly2sym (P,p);
            return
        end
end

%% (connected) graph with node of n-1 degree
%
if (max(Es) == n-1)  % may be optimized for multiple nodes
    cntMAXDEG = cntMAXDEG + 1;
    max_mask = Es==(n-1);
    disp(['[INFO] found ',int2str(nnz(max_mask)),' nodes of degree n-1!'])

    if true %nnz(max_mask) == 1 % single node
        maxdeg = find(Es==(n-1),1);
        mask = [1:maxdeg-1 maxdeg+1:n];
     
        %m = V(maxdeg)
        %Wpol  = sym2poly_arr(W);
        %VWpol = sym2poly_arr(VW); can't do this because of unpredictable result length!
        %VWpol = Wpol2VWpol(rel,Wpol);
        %VW(maxdeg)  = W(maxdeg); % VW will be used only if maxdeg-th node works;
        VWpol(maxdeg,:) = [0 Wpol(maxdeg,:)];
        
        if  rel(maxdeg) %m == 1
            P = ECPN_full_pol(VWpol);
        else            %m == p
            P = [ECPN_full_pol(VWpol) 0]; %P = p * ECPN_full(VW);
            tmp = ECPN_pol(rel(mask),E(mask,mask),Wpol(mask,:));
            tmp = conv([-1,1],tmp);
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
    P = poly2sym(P,p);
    return
end

%% Debug -- checking for hanging trees(or nodes or chains)
%
%also should find trees really good
if ~isempty(find(Es==1,1))
    cntHNODES = cntHNODES + 1;
%    disp('Found hnodes!');
    if (q == n-1)
%    cntTREE = cntTREE + 1;
        disp(['Found tree = ' int2str(n)])
    end
    %P = ECPN_hnodes(V,E,W);
    P = ECPN_hnodes_pol(rel,E,Wpol);
    P = poly2sym(P,p);
    return;

end

%% Reduction of remaining chains (if any)
%

%detect chains
ind_ch = Es==2;
    
%have to temporarily disable this optimization
%because of conflicts with reliable nodes reducing
%to enable ASAP
if nnz(ind_ch) > 0 %have 2 or more nodes with deg = 2
    E_ch = E;
    E_ch(~ind_ch,:) = 0;
    E_ch(:,~ind_ch) = 0;
    [CompNum comps] = graphconncomp(E_ch,'Directed',false);
    if ~(CompNum < n)
        disp('0 chain(s) found!');
    else % have chain
        c_len = zeros(0,CompNum);
        for i = 1:CompNum
            c_len(i) = nnz(comps == i); 
        end
        c_len(c_len == 1) = 0;
    
        [c_max_len c_max] = max(c_len);

%        disp([int2str(nnz(c_len)), ' chain(s) found (max length = ', int2str(c_max_len), '): ', int2str(nonzeros(c_len)')]);
        
        c = find(comps == c_max);
        if length(c) > 2
            Es_ch = sum(E_ch);
            c = c(graphtraverse(E(c,c),find(Es_ch(c)==1,1),'Directed',false)); %traversing chain
        end
        c0   = find(E(c(1),:)   & ~E_ch(c(1),:)  );
        cend = find(E(c(end),:) & ~E_ch(c(end),:));
        index = 1:n;
        
        if c0 == cend 
            disp('OLOLO! c0 = cend!');
        end
        
        index([c0 c]) = [];
        index = [c0 c index];
        cend = find(index == cend);
        
        %reordering nodes, maybe need to avoid this
        Wpol=sym2poly_arr(W);
        Wpol  =  Wpol(index,:); 
        %V  =  V(index);
        %VW = VW(index);
        rel = rel(index);
        E  =  E(index,index);
        
        %Es = Es(index);
        %now we have biggest chain on nodes 2:c_max_len+1
        %nodes 1:2 and (c_max_len+1):(c_max_len+2) are adjacent
        
        %need to optimize chain reductions
        
        Qpol = rel2Qpol(rel);
        
        %mult = prod(V(2:c_max_len+1));
        mult = [1 zeros(1,nnz(~rel(2:c_max_len+1)))];
%        mult = poly2sym(mult,p);
        P = 0;
        mult_pol = [1];
        
        cum_nonrel = [0,cumsum(~rel(2:end))]; % look into making array shorter
        for i=2:c_max_len+1
        
            ind_q = [1:i-1 i+1:n];            
            
            %P = P + mult*Q(i) * ECPN(V(ind_q),E(ind_q,ind_q),W(ind_q));
                       
            tmp = ECPN_pol(rel(ind_q),E(ind_q,ind_q),Wpol(ind_q,:));
            tmp = conv( conv(mult_pol,Qpol(i,:)) , tmp );
            P = poly_add(P,tmp);
            
            mult_pol = [1 zeros(1,cum_nonrel(i))];
            rel(i) = 1;
        end
        
        %P = P + mult*ECPN_full(W(2:c_max_len+1));
         
        tmp = conv( mult_pol , ECPN_full_pol(Wpol(2:c_max_len+1,:)) );
        P = poly_add(P,tmp);
        
        Wpol(2,:) = sum(Wpol(2:c_max_len+1,:),1); %!!! CAREFUL !!!

%        if cend ~= 2
            E(2,cend) = 1;
            E(cend,2) = 1;
%        else 
%            disp('OLOLO! cend = 2');
%        end
        ind_p = [1:2 c_max_len+2:n];
        %gcc = graphconncomp(E(ind_p,ind_p),'Directed',false);
        %if gcc > 1
        %    disp(['conncomp = ', int2str()]);
        %end
        V = rel2V(rel);
        for i = n:-1:1
            W(i)  = poly2sym(Wpol(i,:),p);
        end
        %P = poly2sym(P,p);
        %P = P + mult*ECPN_C(V(ind_p),E(ind_p,ind_p),W(ind_p));
        tmp = ECPN_C(V(ind_p),E(ind_p,ind_p),W(ind_p));
        tmp = sym2poly(tmp);
        P = poly_add(P, conv(mult, tmp) );
        P = poly2sym(P,p);
        return
    end
end

%% Fully reliable connected graph
%
%think about moving to ECPN
assert(nnz(~rel)>0,'[ERROR] Assertion failed -- found fully reliable connected graph in ECPN_C just before branching!');
if nnz(~rel)==0 %V == ones(1,n) %sum(V) == n
    cntRELIABLE = cntRELIABLE + 1;
    P = ECPN_full_pol(Wpol); %using W instead of WV
    P = poly2sym(P,p);
    return
end

%%
%
cntBRANCHING = cntBRANCHING + 1;

%%% [V_MaxDegree R] = max(EIs);
%%% R = R(1);
%%% C = find(EI(R,:),1); % but maybe it's better to find another high-degree node

%Vnonrel = find((V ~=0)&(V ~= 1));
%red = Vnonrel(1);


%% Sorting nodes by degree; do not need this ATM
%

[~, index] = sort(Es,'descend');
%[Es, index] = sort(Es,'descend');
%index = colperm(E);
%Es = Es(index);

rel = rel(index);
E = E(index,index);
%W = W(index);
Wpol = Wpol(index,:);
%VW = VW(index);

%ascending order slows up things!
%needs testing

%% choosing node for reduction
%

%red = find((V~=0)&(V~=1),1);
%red = find((V~=1),1); % no zero nodes anyway ATM
red = find(~rel,1); % no zero nodes anyway ATM

%red = 1; % very-very QnD hack; needs testing
%if V(red)==1
%    disp('CRITICAL ASSUMPTION FAILURE!!!!!!!!!!!!!!');
%    red = find((V~=1),1); % no zero nodes anyway ATM
%end

%[~, index] = sort(Es,'descend'); %avoiding full sotring
%red = index(find((V(index)~=1),1));
%red = index(1);
%disp(int2str(red));

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

%% first subgraph
ind = true(1,n);
ind(red) = 0;

P = conv(ECPN_pol(rel(ind),E(ind,ind),Wpol(ind,:)),[-1, 1]);
%P = P + (1-m) * ECPN(Vred, Ered, Wred);


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

%% let's contract E without deleting red node
E = ContractNodes(E,red);

%% second subgraph
P = poly2sym(P,p);
V = rel2V(rel);
for i = n:-1:1
    W(i)  = poly2sym(Wpol(i,:),p);
end

%P = P + p * ECPN_C(V(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),W(~nodes_rel_ind));
tmp = ECPN_C(V(~nodes_rel_ind),E(~nodes_rel_ind,~nodes_rel_ind),W(~nodes_rel_ind));
tmp = sym2poly(tmp);
P = poly_add(P, [tmp 0]);
P = poly2sym(P,p);

%removing nodes_rel;
%V(nodes_rel) = [];
%E(nodes_rel,:) = []; E(:,nodes_rel) = [];
%W(nodes_rel) = [];
%P = P + p * ECPN_C(V,E,W);
%P = (1-m) * ECPN(V(ind), E(ind,ind), W(ind)) + m * ECPN_C(V,E,W);

end