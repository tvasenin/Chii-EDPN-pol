function P = ECPN_pol(rel,E,Wpol, cut_idx)
%ECPN Calculates ECPN for given G(V,E,W)
%   Detailed explanation goes here

global cnt

%% Precalc and delete isolated nodes or nodes without weight
%
Es = sum(E);
Vzero = find(Es==0);

if ~isempty(Vzero)
    rel(Vzero)   = [];
    if isempty(rel) %no edges
        cnt.TOTAL = cnt.TOTAL + 1;
        cnt.NOEDGES = cnt.NOEDGES + 1;
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
[V_comp, c_lens] = components_mex(E); % shortcut, need to place MEX-file to the MATLAB path
CompNum = length(c_lens);

if CompNum > 1
    cnt.TOTAL = cnt.TOTAL + 1;
    cnt.MULTICOMP = cnt.MULTICOMP + 1;
    P = 0;
    for i = CompNum:-1:1
        mask = find(V_comp == i);
        %P = P + ECPN_C(V(mask),E(mask,mask),W(mask));
        if nargin < 4
            tmp = ECPN_C_pol(rel(mask),E(mask,mask),Wpol(mask,:)); % do not calculate cut_idx if it's not necessary
        else
            tmp = ECPN_C_pol(rel(mask),E(mask,mask),Wpol(mask,:), cut_idx(mask)); %rely on correct cut_idx from argin
        end
        P = poly_add(P,tmp);
        %ECPN_conn(i) = ECPN(rel(mask),E(mask,mask),W(mask)); %bad typecast
    end
    %P = sum(ECPN_conn);
else
    if nargin < 4
        P = ECPN_C_pol(rel,E,Wpol); % do not calculate cut_idx if it's not necessary
    else
        P = ECPN_C_pol(rel,E,Wpol, cut_idx); % rely on correct cut_idx from argin
    end
end

end