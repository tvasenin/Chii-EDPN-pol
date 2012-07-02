function P = ECPN_pol(rel,E,Wpol)
%ECPN Calculates ECPN for given G(V,E,W)
%   Detailed explanation goes here

global cntNOEDGES   cntMULTICOMP    cntTOTAL

%% Precalc and delete isolated nodes or nodes without weight
%
Es = sum(E);
Vzero = find(Es==0);

if ~isempty(Vzero)
    rel(Vzero)   = [];
    if isempty(rel) %no edges
        cntTOTAL = cntTOTAL + 1;
        cntNOEDGES = cntNOEDGES + 1;
        P = 0;
        return;
    end
    E(Vzero,:) = [];
    E(:,Vzero) = [];
    Wpol(Vzero,:)   = [];
end

%% From now on we do not have isolated nodes or empty graph.
% Will be careful not to get them from now on.
% TODO: insert assert

%% Connection components count
%[CompNum V_comp] = graphconncomp(E,'Directed',false);
[CompNum V_comp] = graphalgs('wcc',0,false,E); % shortcut, need to place graphalgs MEX-file to the MATLAB path

if CompNum > 1
    cntTOTAL = cntTOTAL + 1;
    cntMULTICOMP = cntMULTICOMP + 1;
    P = 0;
    for i = CompNum:-1:1
        mask = find(V_comp == i);
        %P = P + ECPN_C(V(mask),E(mask,mask),W(mask));
        tmp = ECPN_C_pol(rel(mask),E(mask,mask),Wpol(mask,:));
        P = poly_add(P,tmp);
        %ECPN_conn(i) = ECPN(rel(mask),E(mask,mask),W(mask)); %bad typecast
    end
    %P = sum(ECPN_conn);
else
P = ECPN_C_pol(rel,E,Wpol);
end

end