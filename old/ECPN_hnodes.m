function P = ECPN_hnodes(V,E,W)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here

rel = (V == 1);
%VWpol = sym2poly_arr(VW);

Es = sum(E);
VW = V.*W;
n = numel(rel);
q = sum(Es)/2;
P = 0;

hnodes = find(Es==1); %maybe it's possible to rewrite via hnodes = Es==1;
ind = true(1,n);
    
%disp('===================================')
while ~isempty(hnodes)
    if q == 1 %have only 2 nodes
%            P = P + ECPN_C_numel2(V(hnodes).*W(hnodes));
        P = P + ECPN_C_numel2(VW(hnodes)); % now we recalc VW.
        q = 0;
        break
    end
%   if q == 3;
%       P = ECPN_C_numel3(V(hnodes),E,VW);
%       q = 0;
%       break;
%   end

    %switch length(hnodes)
    [i1, hneis] = find(E(hnodes,:));
    [~, i2] = sort(i1); %to optimize!
    hneis = hneis(i2);

    [hneis, ind2, ~] = unique(hneis); % to optimize
    hnodes = hnodes(ind2);
    
 %      if Es(hnei) == 2
 %          disp('WOW! Hanging chain has been found!');
 %      end
    P = P + VW(hnodes)*VW(hneis)'; %shouldn't use VW for hnode!!
    W(hneis) = W(hneis) + VW(hnodes);
    %Es(hneis)  = Es(hneis)  - 1;
    %Es(hnodes) = Es(hnodes) - 1; %just to be safe
        
%   q = q - length(hnodes);
    E(hnodes,:) = 0;
    E(:,hnodes) = 0;
    ind(hnodes) = 0;
    Es = sum(E);
    q = sum(Es)/2;

    VW = V.*W; % recalc VW; may be faster than sym.subsref
    hnodes = find(Es==1);
end
        
if q > 0 %still have some edges; rewrite with q ~= 0
    P = P + ECPN_C(V(ind),E(ind,ind),W(ind));
end

%syms p
%P = poly2sym(P,p);

end