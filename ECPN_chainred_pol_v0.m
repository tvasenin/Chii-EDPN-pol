function P = ECPN_chainred_pol_v0(rel,E,Wpol)
%ECPN_chainred_pol Calculates ECPN for a graph with chain inside
%   Detailed explanation goes here

%% assert: no hnodes, not a cycle

n = numel(rel);
Es = sum(E);
P = zeros(1,2*length(Wpol(1,:))+1);
ind_ch = (Es == 2);


E_ch = E;
E_ch(~ind_ch,:) = 0;
E_ch(:,~ind_ch) = 0;
%[CompNum comps] = graphconncomp(E_ch,'Directed',false);
[CompNum comps] = graphalgs('wcc',0,false,E_ch); % shortcut, need to place graphalgs MEX-file to the MATLAB path
%assert((CompNum < n) == (nnz(E(ind_ch,ind_ch)) > 0), 'OLOLO');

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
    %            c = c(graphtraverse(E(c,c),find(Es_ch(c)==1,1),'Directed',false)); %traversing chain
    c = c(graphalgs('dfs',0,false,E(c,c),find(Es_ch(c)==1,1),inf));% shortcut, need to place graphalgs MEX-file to the MATLAB path
end
c0   = find(E(c(1),:)   & ~E_ch(c(1),:)  );
cend = find(E(c(end),:) & ~E_ch(c(end),:));
index = 1:n;

%if c0 == cend
%    disp('OLOLO! c0 = cend!');
%end

index([c0 c]) = [];
index = [c0 c index];
cend = find(index == cend);

%reordering nodes, maybe need to avoid this
Wpol  =  Wpol(index,:);
%V  =  V(index);
%VW = VW(index);
rel = rel(index);
E  =  E(index,index);

%Es = Es(index);
%now we have biggest chain on nodes 2:c_max_len+1
%nodes 1:2 and (c_max_len+1):(c_max_len+2) are adjacent

%need to optimize chain reductions

%mult = prod(V(2:c_max_len+1));
%mult = [1 zeros(1,nnz(~rel(2:c_max_len+1)))];
%mult_len = nnz(~rel(2:c_max_len+1));
cum_nonrel = [0,cumsum(~rel(2:end))]; % look into making array shorter

P = 0;
%         mult_pol = [1];
%[chainbridge comps2([1 3:n])] = graphconncomp(E([1 3:n],[1 3:n]),'Directed',false);
%[chainbridge comps_br([1 3:n])] = graphalgs('wcc',0,false,E([1 3:n],[1 3:n])); % shortcut, need to place graphalgs MEX-file to the MATLAB path
%comps_br(2) = 1; % important!
[chainbridge comps_br([1 c_max_len+2:n])] = graphalgs('wcc',0,false,E([1 c_max_len+2:n],[1 c_max_len+2:n])); % shortcut, need to place graphalgs MEX-file to the MATLAB path
assert(chainbridge < 3, '[ERROR] Assertion failed: Unbelieveable! First node of the chain breaks graph into >2 conn comps!');
assert(comps_br(1) == 1,'[ERROR] Left part of grapg has index 2!');

ind1 = (comps_br == 1);
ind2 = (comps_br == 2);
%        cnext = [2:c_max_len+1 cend];

if ~rel(2) % unfolding loop
    ind_q = true(1,n);
    ind_q(2)= false;
%    if chainbridge == 2
%        ind_q_l = ind1;
%        ind_q_r = ind2;
%        ind_q_r(3:c_max_len+1) = true;
%        
%        %ind_q_l = ind1 & ind_q;
%        %ind_q_r = ind2 & ind_q;
%        tmp =               ECPN_C_pol(rel(ind_q_l),E(ind_q_l,ind_q_l),Wpol(ind_q_l,:));
%        tmp = poly_add(tmp, ECPN_C_pol(rel(ind_q_r),E(ind_q_r,ind_q_r),Wpol(ind_q_r,:)));
%    else
%        tmp = ECPN_C_pol(rel(ind_q),E(ind_q,ind_q),Wpol(ind_q,:));
        tmp = ECPN_pol(rel(ind_q),E(ind_q,ind_q),Wpol(ind_q,:));
%    end
    tmp = -[tmp 0] + [0 tmp];
    P = poly_add(P,[tmp  zeros(1,cum_nonrel(1))]);
    rel(2) = 1;
end

Wpol2 = Wpol(2,:);
tmp_c = zeros(1,2*length(Wpol(1,:))-1); % for elimination of poly_add
%now we can reduce resulting reliable node in chain to node 2 in each iteration
for i=3:c_max_len+1 % add one to c_max_len because of the starting node!
    if ~rel(i)
        %ind_q = [1:i-1 i+1:n];
        ind_q = true(1,n);
        ind_q(i)= false;
        %                ind_q(3:i)= false;
        %P = P + mult*Q(i) * ECPN(V(ind_q),E(ind_q,ind_q),W(ind_q));
        
%        if chainbridge == 2
%            %ind_q_l = ind1 & ind_q;
%            %ind_q_r = ind2 & ind_q;
%            ind_q_l = ind1;
%            ind_q_l(2:i-1) = true;
%            ind_q_r = ind2;
%            ind_q_r(i+1:c_max_len+1) = true;
%            tmp =               ECPN_C_pol(rel(ind_q_l),E(ind_q_l,ind_q_l),Wpol(ind_q_l,:));
%            tmp = poly_add(tmp, ECPN_C_pol(rel(ind_q_r),E(ind_q_r,ind_q_r),Wpol(ind_q_r,:)));
%        else
%            tmp = ECPN_C_pol(rel(ind_q),E(ind_q,ind_q),Wpol(ind_q,:));
            tmp = ECPN_pol(rel(ind_q),E(ind_q,ind_q),Wpol(ind_q,:));
%        end
        %                tmp = poly_add(tmp, tmp_c);
        tmp = -[tmp 0] + [0 tmp];
        P = poly_add(P,[tmp  zeros(1,cum_nonrel(i-1))]);
        rel(i) = 1;
    end
    %            tmp_c = tmp_c + conv2(Wpol(2,:),Wpol(i,:));
%    tmp_c = tmp_c + conv2(Wpol2,Wpol(i,:));
    %P = poly_add(P,[conv2(Wpol(2,:),Wpol(i,:)) zeros(1,cum_nonrel(i-1))]);
%    Wpol2 = Wpol2 + Wpol(i,:);
    %            Wpol(2,:) = Wpol(2,:) + Wpol(i,:);
    %            E(2,cnext(i)) = 1;
    %            E(cnext(i),2) = 1;
end

%P = P + mult*ECPN_full(W(2:c_max_len+1));
        tmp = ECPN_full_pol(Wpol(2:c_max_len+1,:));
        P = poly_add(P,[tmp  zeros(1,cum_nonrel(c_max_len+1))]);
%P = poly_add(P,[tmp_c  zeros(1,cum_nonrel(c_max_len+1))]);
%P = poly_add(P,tmp_c);

        Wpol(2,:) = sum(Wpol(2:c_max_len+1,:),1); %!!! CAREFUL !!!
%Wpol(2,:) = Wpol2;

%        assert (cend ~= 2, 'Assertion failed: unbelieveable, cend ~=2 !');
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

%P = P + mult*ECPN_C(V(ind_p),E(ind_p,ind_p),W(ind_p));
tmp = ECPN_C_pol(rel(ind_p),E(ind_p,ind_p),Wpol(ind_p,:));
P = poly_add(P,[tmp  zeros(1,cum_nonrel(c_max_len+1))]);


%P = P(find(P,1):end); % trim leading zeros
%P = polytrim_fast(P); % trim leading zeros

end