function P = ECPN_C_numel4(V,E,W)
%ECPN_C_numel4 Calculates ECPN for connected graph, where n = 4;
%   Detailed explanation goes here

Es = sum(E);
[sgn,index]= sort(Es,'descend')

q = sum(Es)/2; %number of edges

switch sgn
    

if (q == 6) % [3-3-3-3] special case !
    VWpol = Wpol2VWpol(rel,Wpol);
    P = ECPN_full_pol(VWpol);
    P = poly2sym(P,p);
    return
end

rel = (V == 1);
Wpol = sym2poly_arr(W);
            
syms p
    
switch max(Es)
    case 2 %[2-2-1-1] | [2-2-2-2]
        switch q
            case 3 %[2-2-1-1] 
                VWpol = Wpol2VWpol(rel,Wpol);
                P = ECPN_chain_(V,E,V.*W);
                P = sym2poly(O)
            case 4 %[2-2-2-2]
                P = ECPN_cycle_pol(rel,E,Wpol);
                P = poly2sym(P,p);
            otherwise 
                error('[ERROR] Unknown error in ECPN_C_numel4!');
        end
    case 3 %[3-1-1-1] | [3-2-2-1] | [3-3-2-2]   (no [3-3-3-3] )
        maxdeg = find(Es==3,1);
        rel_maxdeg = rel(maxdeg);
        rel(maxdeg) = 1;
        VWpol = Wpol2VWpol(rel,Wpol);
        P = ECPN_full_pol(VWpol);
        if rel_maxdeg
            return 
        end
        P = conv([1 0], P);
        switch q
            case 3 %[3-1-1-1]
                %P = m * feval(symengine, 'ECPN_full_mu',V.*W);
                tmp = 0;
            case 4 %[3-2-2-1]
                %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel2(V(ind_2).*W(ind_2));
                mask = Es==2;
                tmp = ECPN_C_numel2_pol(VWpol(mask,:));
            case 5 %[3-3-2-2]
                mask = [1:maxdeg-1 maxdeg+1:4];
                %P = m * feval(symengine, 'ECPN_full_mu',V.*W) + (1-m) * ECPN_C_numel3(V(mask),E(mask,mask),W(mask));
                tmp = ECPN_C_numel3_pol(rel(mask),E(mask,mask),Wpol(mask,:));
            case 6
                %had earlier termination
            otherwise
                error('[ERROR] Unknown error in ECPN_C_numel4!');
        end
        P = poly_add(P,conv([-1 1],tmp));
    otherwise
        error('[ERROR] Calling ECPN_C_numel4 for non-connected graph!');
end
    







switch q
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
        
        P = ECPN_full_pol(VWpol);
        if rel_maxdeg
            
        else
            P = poly_add(conv([1 0], P), conv([-1 1], ECPN_C_numel3_pol(rel(mask),E(mask,mask),Wpol(mask,:))) );
        end
            P = poly2sym(P,p);
    case 6
        VWpol = Wpol2VWpol(rel,Wpol);
        P = ECPN_full_pol(VWpol);
        P = poly2sym(P,p);
end
end