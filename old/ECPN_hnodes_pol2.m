function P = ECPN_hnodes_pol2(rel,E,Wpol)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here


n = numel(rel);
%syms p
VWpol = Wpol2VWpol(rel, Wpol); % optimize the routine to not use VW at all!

%Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^

Es = sum(E);
%q = sum(Es)/2;
P = 0;

ind = true(1,n);
hnodes = (Es == 1);
q = sum(Es)/2;
%disp('===================================')
while nnz(hnodes)
    switch q
        case 0
            return
        case 1  %have only 2 nodes
            %P = P + ECPN_C_numel2(VW(hnodes)); % now we recalc VW.
            P = poly_add(P,ECPN_C_numel2_pol(VWpol(hnodes,:))); % now we recalc VW.
            %VWpol_hnodes = Wpol2VWpol(rel(hnodes), Wpol(hnodes,:));
            %P = poly_add(P,ECPN_C_numel2_pol(VWpol_hnodes)); % now we recalc VW.
%            q = 0;
%            break
            return
%        case 3
%            P = ECPN_C_numel3 (V(hnodes),E,VW);
%            q = 0;
%            break;
%            Es(:) = 0;
%            hnodes = [];            %break
        otherwise
            for i = 1:length(hnodes)
                if 
            end
            %switch length(hnodes)
            [i1, hneis] = find(E(hn,:));
            [~, i2] = sort(i1); %to optimize!
            hneis = hneis(i2);

            %[hneis, ind2, ~] = unique(hneis); % to optimize
            %[hneis, ind2] = unique_fast(hneis); % to optimize
            
            hnodes = hnodes(ind2);
    
%           if Es(hnei) == 2
%               disp('WOW! Hanging chain has been found!');
%           end
           %P = P + VW(hnodes)*VW(hneis)'; %shouldn't use VW for hnode!!
           %tmp = zeros(1,length(hnodes));
           
           tmp = conv2( VWpol(hnodes(1),:) , VWpol(hneis(1),:) );
           for i = 2:length(hnodes)
               tmp = tmp + conv2( VWpol(hnodes(i),:) , VWpol(hneis(i),:) );
           end
           P = poly_add(P,tmp);
           
           %W(hneis) = W(hneis) + VW(hnodes);
           %Wpol = padarray (Wpol,[0,1],'pre'); % now dim(Wpol) should be equal to dim(VWpol)
           %WVpol_hnodes = VWpol(hnodes,:);
           %if 
           Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
    
           Wpol(hneis,:) = Wpol(hneis,:) + VWpol(hnodes,:);
           %VW = V.*W; % recalc VW; may be faster than sym.subsref
           VWpol = Wpol2VWpol(rel, Wpol); %now dims are not equal again
           
           %VWpol(hneis(rel(hneis))) = 
    
           E(hnodes,:) = 0;
           E(:,hnodes) = 0;
           ind(hnodes) = 0;
           %Es(hneis)  = Es(hneis)  - 1;
           %Es(hnodes) = Es(hnodes) - 1; %just to be safe
           %q = q - length(hnodes);
           Es = sum(E);
           %q = sum(Es)/2;

           hnodes = find(Es==1);
    end
end
        
%V = rel2V(rel);
%for i = n:-1:1
%    W(i)  = poly2sym(Wpol(i,:),p);
%end

tmp = ECPN_C_pol(rel(ind),E(ind,ind),Wpol(ind,:));
P = poly_add(P,tmp);

end