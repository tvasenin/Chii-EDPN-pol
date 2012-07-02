function P = ECPN_chain(V,E,VW)
%ECPN_chain Calculates ECPN for a (single) chain
%   Detailed explanation goes here

rel = (V == 1);
VWpol = sym2poly_arr(VW);

switch numel(V)
    case 1
        disp('[WARNING] Calling ECPN_chain with n = 1');
        P = 0;
    case 2
        disp('[WARNING] Calling ECPN_chain with n = 2');
        P = ECPN_C_numel2_pol(VWpol);
    case 3
        disp('[WARNING] Calling ECPN_chain with n = 3');
        P = ECPN_C_numel3_pol(rel,E,VWpol);
    otherwise
        s = find(sum(E)==1,1); %finding starting node
        index =  graphtraverse(E,s,'Directed',false); %traversing chain
%        P = ECPN_ordered_chain(V(index),W(index),V(index).*W(index)); % resorting (don't need E and VW anymore)
%        P = ECPN_ordered_chain_mu(V(index),W(index),VW(index)); % resorting (don't need E and VW anymore)
        P = ECPN_ordered_chain_pol(rel(index),VWpol(index,:)); % resorting (don't need E and VW anymore)
        
%       P = feval(symengine, 'ECPN_ordered_chain_mu',V(index),VW(index)); % resorting (don't need E and VW anymore)
end

syms p
P = poly2sym(P,p);

end