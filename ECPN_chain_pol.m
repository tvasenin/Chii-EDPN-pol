function P = ECPN_chain_pol(rel,E,VWpol)
%ECPN_chain Calculates ECPN for a (single) chain
%   Detailed explanation goes here

switch numel(rel)
    case 1
        disp('[WARNING] Calling ECPN_chain with n = 1');
        P = 0;
    case 2
        disp('[WARNING] Calling ECPN_chain with n = 2');
        P = ECPN_C_numel2_pol(VWpol); % cant use v2 because it needs rel and Wpol
    case 3
        disp('[WARNING] Calling ECPN_chain with n = 3');
        P = ECPN_C_numel3_pol(rel,E,VWpol);
    otherwise
        index = chaintraverse_fast(E);
        P = ECPN_ordered_chain_pol(rel(index),VWpol(index,:)); % resorting (don't need E and VW anymore)
end

end