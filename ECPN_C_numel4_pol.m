function P = ECPN_C_numel4_pol(rel,E,Wpol)
%ECPN_C_numel4 Calculates ECPN for connected graph, where n = 4;
%   Detailed explanation goes here

Es = sum(E);
q = sum(Es)/2; %number of edges

switch q
    case 3
        if max(Es)==3 %[3-1-1-1]
            maxdeg = find(Es==3,1);
            rel_maxdeg = rel(maxdeg);
            rel(maxdeg) = 1;
            VWpol = Wpol2VWpol(rel,Wpol);
            
            %P = m * feval(symengine, 'ECPN_full_mu',V.*W);
            P = ECPN_full_pol(VWpol);
            if ~rel_maxdeg
                P = [P 0];
            end
        else          %[2-2-1-1]
            VWpol = Wpol2VWpol(rel,Wpol);
            P = ECPN_chain_pol(rel,E,VWpol);
        end
    case 4
        if max(Es)==3 %[3-2-2-1]
            maxdeg = find(Es==3,1);
            mask  = Es==2;
            
            rel_maxdeg = rel(maxdeg);
            rel(maxdeg) = 1;
            VWpol = Wpol2VWpol(rel,Wpol);

            %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel2(V(ind_2).*W(ind_2));
            P = ECPN_full_pol(VWpol);
            if ~rel_maxdeg
                %P = poly_add( [P 0], conv2([-1 1], ECPN_C_numel2_pol(VWpol(mask,:))));
                P = poly_add( [P 0], conv2([-1 1], ECPN_C_numel2_pol_v2(rel(mask),Wpol(mask,:))));
%                tmp = ECPN_C_numel2_pol(VWpol(mask,:));
%                P = poly_add( [P 0], -[tmp 0] + [0 tmp]);
            end
        else          %[2-2-2-2]
            P = ECPN_cycle_pol_v2(rel,E,Wpol);
        end
    case 5            %[3-3-2-2]
        maxdeg = find(Es==3,1);
        mask = [1:maxdeg-1 maxdeg+1:4];
        
        rel_maxdeg = rel(maxdeg);
        rel(maxdeg) = 1; % it will be used only if maxdeg-th node works;
        VWpol = Wpol2VWpol(rel,Wpol);
        
        %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel3(V(mask),E(mask,mask),V(mask).*W(mask));
        P = ECPN_full_pol(VWpol);
        if ~rel_maxdeg
            P = poly_add( [P 0], conv2([-1 1], ECPN_C_numel3_pol(rel(mask),E(mask,mask),VWpol(mask,:))) );
%            tmp = ECPN_C_numel3_pol(rel(mask),E(mask,mask),VWpol(mask,:));
%            P = poly_add( [P 0], [-tmp 0] + [0 tmp]);
        end
    case 6            %[3-3-3-3]
        VWpol = Wpol2VWpol(rel,Wpol);
        P = ECPN_full_pol(VWpol);
    otherwise
        disp('[ERROR] Calling ECPN_C_numel4 with q = ',q);
        error('[HALT] Fatal error, termination in ECPN_C_numel4');
end

end