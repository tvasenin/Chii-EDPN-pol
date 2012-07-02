function P = ECPN_chain_pol(rel,E,VWpol)
%ECPN_chain Calculates ECPN for a (single) chain
%   Detailed explanation goes here

switch numel(rel)
    case 1
        disp('[WARNING] Calling ECPN_chain with n = 1');
        P = 0;
    case 2
        disp('[WARNING] Calling ECPN_chain with n = 2');
        P = ECPN_C_numel2_pol(VW);
    case 3
        disp('[WARNING] Calling ECPN_chain with n = 3');
        P = ECPN_C_numel3_pol(rel,E,VWpol);
    otherwise
        %s = find(sum(E)==1,1); %finding starting node
        %index =  graphtraverse(E,s,'Directed',false); %traversing chain
        %index = graphalgs('dfs',0,false,E,find(sum(E)==1,1),inf);% shortcut, need to place graphalgs MEX-file to the MATLAB path
        index = chaintraverse_fast(E);
%       
%        P = ECPN_ordered_chain(V(index),W(index),V(index).*W(index)); % resorting (don't need E and VW anymore)
%        P = ECPN_ordered_chain_mu(V(index),W(index),VW(index)); % resorting (don't need E and VW anymore)
        P = ECPN_ordered_chain_pol(rel(index),VWpol(index,:)); % resorting (don't need E and VW anymore)
        
%       P = feval(symengine, 'ECPN_ordered_chain_mu',V(index),VW(index)); % resorting (don't need E and VW anymore)
end

end