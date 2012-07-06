function P = ECPN_chainred_pol_v3(rel,E,Wpol)
%ECPN_chainred_pol Calculates ECPN for a graph with chain inside
%   Detailed explanation goes here

%% assert: no hnodes, not a cycle

n = length(E);

[ index, c0, cend, c_first, c_last ] = GetIndexedChain(E);

%reordering nodes, maybe need to avoid this
rel = rel(index);
E  =  E(index,index);
Wpol  =  Wpol(index,:);

%if nnz(rel(c_first:c_last)) ~= 0
%    BGView = biograph(triu(E));
%    set(BGView,'ShowArrows','off');%,'LayoutType','equilibrium');
%    view(BGView)
%    rel
%end
%assert(nnz(rel(c_first:c_last)) == 0, '[ERROR] have reliable nodes in the chain!');
%Es = Es(index);
%now we have biggest chain on nodes c_first:c_last
%nodes c0:c_first and (c_last):cend are adjacent

%need to optimize chain reductions

%mult = prod(V(c_first:c_last));
%mult = [1 zeros(1,nnz(~rel(c_first:c_last)))];
%cum_nonrel = [0,cumsum(~rel(c_first:end))]; % look into making array shorter
cum_nonrel(c_first:n)  = cumsum(~rel(c_first:n)); % look into making array shorter

[chainbridge comps_br([c0 c_first c_last+1:n])] = graphalgs('wcc',0,false,E([c0 c_first c_last+1:n],[c0 c_first c_last+1:n])); % shortcut, need to place graphalgs MEX-file to the MATLAB path
comps_br(c_first) = 0; % important!
assert(chainbridge < 3, '[ERROR] Assertion failed: Unbelieveable! First node of the chain breaks graph into >2 conn comps!');
assert(comps_br(c0) == 1,'[ERROR] Left part of graph has index 2!');

%cnext = [2:c_last cend];

% calculating chain 
%cum_r = get_cumsum_ord_chain_right(rel(3:c_last),Wpol(3:c_last,:));


if chainbridge == 2 % c0 ~= cend 
    ind1 = find(comps_br == 1);
    ind2 = find(comps_br == 2);

    if rel(c_first) % unfolding loop
        P = 0;
    else
        ind_q_l = ind1;
        ind_q_r = [c_last ind2];

        tmp =               ECPN_C_pol(rel(ind_q_l),E(ind_q_l,ind_q_l),Wpol(ind_q_l,:));
        
        tmp_q_r = ECPN_ordered_chain_pol(rel(c_first+1:c_last), Wpol2VWpol( rel(c_first+1:c_last), Wpol(c_first+1:c_last,:)));
        tmp = poly_add(tmp,tmp_q_r);
        
        Wlast = get_cum_Wlast(rel(c_first+1:c_last),Wpol(c_first+1:c_last,:));
        k = length(Wlast)-length(Wpol(1,:));
        tmp = poly_add(tmp, ECPN_C_pol(rel(ind_q_r),E(ind_q_r,ind_q_r),[Wlast ; zeros(length(ind_q_r)-1,k) Wpol(ind_q_r(2:end),:)]));
        
        P = -[tmp 0] + [0 tmp];
        rel(c_first) = 1;
    end

    WpolF = Wpol(c_first,:);
    tmp_p = zeros(1,2*length(Wpol(1,:))-1); % for elimination of poly_add
    %now we can reduce resulting reliable node in chain to node 2 in each iteration
    for i=c_first+1:c_last % add one to c_max_len because of the starting node!
        %tmp_p = tmp_p + conv2(Wpol(c_first,:),Wpol(i,:));
        tmp_p = tmp_p + conv2(WpolF,Wpol(i,:));
        %P = poly_add(P,[conv2(Wpol(c_first,:),Wpol(i,:)) zeros(1,cum_nonrel(i-1))]);
        WpolF = WpolF + Wpol(i,:);
%       Wpol(c_first,:) = WpolF;
        %Wpol(c_first,:) = Wpol(c_first,:) + Wpol(i,:);
        if ~rel(i)
            ind_q_l = [ind1(1)  c_first:(i-1)  ind1(2:end)];  
            ind_q_r = [c_last  ind2];

            tmp =               ECPN_C_pol(rel(ind_q_l),E(ind_q_l,ind_q_l),Wpol(ind_q_l,:));
            if i < c_last
                tmp_q_r = ECPN_ordered_chain_pol(rel(i+1:c_last), Wpol2VWpol( rel(i+1:c_last), Wpol(i+1:c_last,:)));
                tmp = poly_add(tmp,tmp_q_r);

                Wlast = get_cum_Wlast(rel(i+1:c_last),Wpol(i+1:c_last,:));
                k = length(Wlast)-length(Wpol(1,:));
                tmp = poly_add(tmp, ECPN_C_pol(rel(ind_q_r),E(ind_q_r,ind_q_r),[Wlast ; zeros(length(ind_q_r)-1,k) Wpol(ind_q_r(2:end),:)]));
            else
                tmp = poly_add(tmp, ECPN_C_pol(rel(ind_q_r(2:end)),E(ind_q_r(2:end),ind_q_r(2:end)), Wpol(ind_q_r(2:end),:)));
            end
            %tmp = poly_add(tmp, tmp_p);
            tmp = -[tmp 0] + [0 tmp];
            P = poly_add(P,[tmp  zeros(1,cum_nonrel(i-1))]);
            rel(i) = 1;
        end
    end

    
    %P = P + mult*ECPN_full(W(c_first:c_last));
    tmp_p = ECPN_ordered_chain_pol(rel(c_first:c_last), Wpol2VWpol( rel(c_first:c_last), Wpol(c_first:c_last,:)));
    P = poly_add(P,[tmp_p  zeros(1,cum_nonrel(c_last))]);
%    P = poly_add(P,tmp_p);
    
    %Wpol(c_first,:) = sum(Wpol(c_first:c_last,:),1); %!!! CAREFUL !!!
    Wpol(c_first,:) = WpolF;
    
    assert (cend ~= c0, 'Assertion failed: unbelieveable, cend==c0 while chainbridge==2 !');
    
    E(c_first,cend) = 1; E(cend,c_first) = 1;
    %need to rewrite the algorhitm to enable next line
%   E(c0,cend) = 1; E(cend,) = 1; % contracting c0 with cend
    ind_p = [c0 c_first c_last+1:n];
    
    %P = P + mult*ECPN_C(V(ind_p),E(ind_p,ind_p),W(ind_p));
    tmp = ECPN_C_pol(rel(ind_p),E(ind_p,ind_p),Wpol(ind_p,:));
    P = poly_add(P,[tmp zeros(1,cum_nonrel(c_last))]);

else % chainbridge = 1
    
    if ~rel(c_first) % unfolding loop
        ind_q = [c0  c_last:n];
        
        Wlast = get_cum_Wlast(rel(c_first+1:c_last),Wpol(c_first+1:c_last,:));
        k = length(Wlast)-length(Wpol(1,:));
        Wcomp = [zeros(length(ind_q),k) Wpol(ind_q,:)];
        Wcomp(2,:) = Wlast;
        tmp = ECPN_C_pol(rel(ind_q),E(ind_q,ind_q),Wcomp);
        tmp_q_r = ECPN_ordered_chain_pol(rel(c_first+1:c_last), Wpol2VWpol( rel(c_first+1:c_last), Wpol(c_first+1:c_last,:)));
        tmp = poly_add(tmp,tmp_q_r);
        
        tmp = -[tmp 0] + [0 tmp];
        P = [tmp  zeros(1,cum_nonrel(1))];
        rel(c_first) = 1;
    else
        P = 0;
    end

    WpolF = Wpol(c_first,:);
    tmp_p = zeros(1,2*length(Wpol(1,:))-1); % for elimination of poly_add
    %now we can reduce resulting reliable node in chain to node 2 in each iteration
    for i=c_first+1:c_last % add one to c_max_len because of the starting node!
        if ~rel(i)

             if i < c_last
                 ind_q = [c0 c_first:i-1 c_last:n];
                 Wlast = get_cum_Wlast(rel(i+1:c_last),Wpol(i+1:c_last,:));
                 k = length(Wlast)-length(Wpol(1,:));
                 Wcomp = [zeros(length(ind_q),k) Wpol(ind_q,:)];
                 Wcomp(i-c_first+2,:) = Wlast;
                 tmp = ECPN_C_pol(rel(ind_q),E(ind_q,ind_q), Wcomp);
                 tmp_q_r = ECPN_ordered_chain_pol(rel(i+1:c_last), Wpol2VWpol( rel(i+1:c_last), Wpol(i+1:c_last,:)));
                 tmp = poly_add(tmp,tmp_q_r);
                 
             else %i == c_last
                 ind_q = [c0 c_first:i-1 c_last+1:n];
                 tmp = ECPN_C_pol(rel(ind_q),E(ind_q,ind_q), Wpol(ind_q,:));
             end

            tmp = -[tmp 0] + [0 tmp];
            P = poly_add(P,[tmp  zeros(1,cum_nonrel(i-1))]);
            rel(i) = 1;
        end
        %tmp_p = tmp_p + conv2(Wpol(2,:),Wpol(i,:));
        tmp_p = tmp_p + conv2(WpolF,Wpol(i,:));
        %P = poly_add(P,[conv2(Wpol(2,:),Wpol(i,:)) zeros(1,cum_nonrel(i-1))]);
        WpolF = WpolF + Wpol(i,:);
        %Wpol(c_first,:) = Wpol(c_first,:) + Wpol(i,:);
        %E(c_first,cnext(i)) = 1;
        %E(cnext(i),c_first) = 1;
    end

    %P = P + mult*ECPN_full(W(c_first:c_last));
    %tmp = ECPN_full_pol(Wpol(c_first:c_last,:));
    %P = poly_add(P,[tmp  zeros(1,cum_nonrel(c_last))]);
    %P = poly_add(P,[tmp_p  zeros(1,cum_nonrel(c_last))]);
    %P = poly_add(P,tmp_p);

    %Wpol(c_first,:) = sum(Wpol(c_first:c_last,:),1); %!!! CAREFUL !!!
    Wpol(c_first,:) = WpolF;

    assert (cend ~= c_first, 'Assertion failed: unbelieveable, cend ~=c_first !');

    if cend == c0 %c0 == cend
        if rel(c0)
   %        tmp_p = poly_add(tmp_p,conv(Wpol(c0,:),Wpol(c_first,:)));
            tmp_p = tmp_p + conv(Wpol(c0,:),Wpol(c_first,:));
        else
            tmp_p = [0 tmp_p] + [conv(Wpol(c0,:),Wpol(c_first,:)) 0];
        end
        Wpol(c0,:) = Wpol(c0,:) + Wpol(c_first,:);
        ind_p = [c0 c_last+1:n];
    else
        E(c_first,cend) = 1; E(cend,c_first) = 1;
        %need to rewrite the algorhitm to enable next line
%        E(c0,cend) = 1; E(cend,c0) = 1; % contracting c0 with cend
        ind_p = [c0 c_first c_last+1:n];
    end
    
    %P = P + mult*ECPN_C(V(ind_p),E(ind_p,ind_p),W(ind_p));
    tmp_p = poly_add(tmp_p,ECPN_C_pol(rel(ind_p),E(ind_p,ind_p),Wpol(ind_p,:)));
    P = poly_add(P,[tmp_p  zeros(1,cum_nonrel(c_last))]);
end

end