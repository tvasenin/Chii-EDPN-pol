function P = ECPN_C_numel4(V,E,W)
%ECPN_C_numel4 Calculates ECPN for connected graph, where n = 4;
%   Detailed explanation goes here

Es = sum(E);
q = sum(Es)/2; %number of edges

rel = (V == 1);
Wpol = sym2poly_arr(W);
            
syms p

%switch max(Es)
%    case 2
%    case 3
%    otherwise
%        disp('[ERROR] Calling ECPN_C_numel4 with q = ',q)
%end
    
switch q
    case 0
        disp('[ERROR] Calling ECPN_C_numel4 with q = 0')
    case 1
        disp('[ERROR] Calling ECPN_C_numel4 with q = 1')
    case 2
        disp('[ERROR] Calling ECPN_C_numel4 with q = 2')
    case 3
        if max(Es)==3 %[3-1-1-1]
            maxdeg = find(Es==3,1);
            m = V(maxdeg);
            V(maxdeg) = 1; 
            P = m * feval(symengine, 'ECPN_full_mu',V.*W);
        else          %[2-2-1-1]
            P = ECPN_chain(V,E,V.*W);
        end
    case 4
        if max(Es)==3 %[3-2-2-1]
            maxdeg = find(Es==3,1);
            %m = V(maxdeg);
	        %V(maxdeg) = 1; 
%            ind_2  = find(Es==2);
            mask  = Es==2;
            
            %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel2(V(ind_2).*W(ind_2));

            rel_maxdeg = rel(maxdeg);
            rel(maxdeg) = 1;
            
            VWpol = Wpol2VWpol(rel,Wpol);
            
            if rel_maxdeg
                P = ECPN_full_pol(VWpol);
            else
                P = conv([1 0], ECPN_full_pol(VWpol));
                P = poly_add(P, conv([-1 1], ECPN_C_numel2_pol(VWpol(mask,:))));
            end
        
            P = poly2sym(P,p);
            
        else          %[2-2-2-2]
            P = ECPN_cycle_pol(rel,E,Wpol);
            P = poly2sym(P,p);
        end
    case 5
        maxdeg = find(Es==3,1);
        mask = [1:maxdeg-1 maxdeg+1:4];
        
        %m = V(maxdeg);
        rel_maxdeg = rel(maxdeg);
        %V(maxdeg) = 1; % it will be used only if maxdeg-th node works;
        rel(maxdeg) = 1; % it will be used only if maxdeg-th node works;
        
        VWpol = Wpol2VWpol(rel,Wpol);
        
        %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel3(V(mask),E(mask,mask),W(mask));
        
        if rel_maxdeg
            P = ECPN_full_pol(VWpol);
        else
            P = conv([1 0], ECPN_full_pol(VWpol));
            P = poly_add(P, conv([-1 1], ECPN_C_numel3_pol(rel(mask),E(mask,mask),Wpol(mask,:))) );
        end
            P = poly2sym(P,p);
    case 6
        VWpol = Wpol2VWpol(rel,Wpol);
        P = ECPN_full_pol(VWpol);
        P = poly2sym(P,p);
end
end