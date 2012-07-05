function P = ECPN_hnodes_pol_v0(rel,E,Wpol)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here

n = numel(rel);
VWpol = Wpol2VWpol(rel, Wpol); % optimize the routine to not use VW at all!

%Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^

Es = sum(E);
P = 0;

hnodes = find(Es==1); %maybe it's possible to rewrite via hnodes = Es==1;
ind = true(1,n);
q = sum(Es)/2;
%disp('===================================')
while ~isempty(hnodes)
    switch q
        case 1  %have only 2 nodes
            %P = P + ECPN_C_numel2(VW(hnodes)); % now we recalc VW.
            P = poly_add(P,ECPN_C_numel2_pol(VWpol(hnodes,:))); % now we recalc VW.
            %VWpol_hnodes = Wpol2VWpol(rel(hnodes), Wpol(hnodes,:));
            %P = poly_add(P,ECPN_C_numel2_pol(VWpol_hnodes)); % now we recalc VW.
%            q = 0;
%            break
            return
        otherwise

            [i1, neis] = find(E(hnodes,:));
            [~, i2] = sort(i1); %to optimize!
            neis = neis(i2);

            %[hneis, ind2, ~] = unique(hneis); % to optimize
            [neis, ind2] = unique_fast(neis); % to optimize
            
            hnodes = hnodes(ind2);
    
%           if Es(hnei) == 2
%               disp('WOW! Hanging chain has been found!');
%           end
           %P = P + VW(hnodes)*VW(hneis)'; %shouldn't use VW for hnode!!

           tmp = conv2( VWpol(hnodes(1),:) , VWpol(neis(1),:) );
           for i = 2:length(hnodes)
               tmp = tmp + conv2( VWpol(hnodes(i),:) , VWpol(neis(i),:) );
           end
           
           P = poly_add(P,tmp);
           
           %W(hneis) = W(hneis) + VW(hnodes);
           %Wpol = padarray (Wpol,[0,1],'pre'); % now dim(Wpol) should be equal to dim(VWpol)
           Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
           Wpol(neis,:) = Wpol(neis,:) + VWpol(hnodes,:);
           %VW = V.*W; % recalc VW; may be faster than sym.subsref
           VWpol = Wpol2VWpol(rel, Wpol); %now dims are not equal again
           
    
           E(hnodes,:) = 0;
           E(:,hnodes) = 0;
           ind(hnodes) = 0;
           %Es(hneis)  = Es(hneis)  - 1;
           %Es(hnodes) = Es(hnodes) - 1; %just to be safe
           %q = q - length(hnodes);
           Es = sum(E);
           q = sum(Es)/2;

           hnodes = find(Es==1);
    end
end


if nnz(Es) > 0 %still have some edges
    tmp = ECPN_C_pol(rel(ind),E(ind,ind),Wpol(ind,:));
    P = poly_add(P,tmp);
end

end