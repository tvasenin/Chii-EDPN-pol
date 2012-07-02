function P = ECPN(V,E,W)
%ECPN Calculates ECPN for given G(V,E,W)
%   Detailed explanation goes here

global cntNOEDGES   cntMULTICOMP    cntTOTAL

%% Precalc and delete isolated nodes or nodes without weight
%

Es = sum(E);
Vzero = find(Es==0);

if ~isempty(Vzero)
    V(Vzero)   = [];
    if isempty(V) %no edges
        cntTOTAL = cntTOTAL + 1;
        cntNOEDGES = cntNOEDGES + 1;
        P = 0;
        return;
    end
    E(Vzero,:) = [];
    E(:,Vzero) = [];
    W(Vzero)   = [];
end

%% From now on we do not have isolated nodes or empty graph.
% Will be careful not to get them from now on.

%% Preparing for 
%
%ind_rel = find(V==1);
%for i = 1:length(ind_rel)
%    cont = find (E(ind_rel(i),:)~=0);
%    E(cont,cont) = 1;
%end
%E = spdiags(sparse(n,1),0,E); %E = E - diag(diag(E));

%% Connection components count
%
[CompNum V_comp] = graphconncomp(E,'Directed',false);

if CompNum > 1
    cntTOTAL = cntTOTAL + 1;
    cntMULTICOMP = cntMULTICOMP + 1;
    P = 0;
    for i = CompNum:-1:1
        mask = find(V_comp == i);
        P = P + ECPN_C(V(mask),E(mask,mask),W(mask));
        %ECPN_conn(i) = ECPN(V(mask),E(mask,mask),W(mask)); %bad typecast
    end
    %P = sum(ECPN_conn);    
    return
end

P = ECPN_C(V,E,W);
return

end