function P = ECPN_hnodes_pol_v2(rel,E,Wpol)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here

n = numel(rel);
Es = sum(E);
P = 0;

hnodes = find(Es==1); %maybe it's possible to rewrite via hnodes = Es==1;
%ind = true(1,n);
q = sum(Es)/2;
%disp('===================================')
while ~isempty(hnodes)
    switch q
        case 1  %have only 2 nodes
            VWpol = Wpol2VWpol(rel, Wpol); % optimize the routine to not use VW at all!
                %P = P + ECPN_C_numel2(VW(hnodes)); % now we recalc VW.
            %P = poly_add(P,ECPN_C_numel2_pol(VWpol(hnodes,:))); % now we recalc VW.
            P = poly_add(P,ECPN_C_numel2_pol(VWpol)); % now we recalc VW.
            return
        otherwise
            %switch length(hnodes)
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
           tmp = zeros(1,2*length(Wpol(1,:))+1);
           
           for i = 1:length(hnodes)
               hnode = hnodes(i);
               nei = neis(i);
               if      rel(hnode) &&  rel(nei) % [1 1]
                   tmp = tmp + [0  0  conv2(Wpol(hnode,:),Wpol(nei,:))];
               elseif ~rel(hnode) && ~rel(nei) % [p p]
                   tmp = tmp + [conv2(Wpol(hnode,:),Wpol(nei,:))  0  0];
               else                            % [0 p] or [p 0]
                   tmp = tmp + [0  conv2(Wpol(hnode,:),Wpol(nei,:))  0];
               end
           end
           
           P = poly_add(P,tmp);
           
           %W(hneis) = W(hneis) + VW(hnodes);
           %Wpol = padarray (Wpol,[0,1],'pre'); % now dim(Wpol) should be equal to dim(VWpol)
           
           %Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
           %Wpol(neis,:) = Wpol(neis,:) + VWpol(hnodes,:);
           
           hnodes_nonrel = hnodes(~rel(hnodes));
           if ~isempty(hnodes_nonrel)
               neis_hnonrel  =   neis(~rel(hnodes));
               hnodes_rel    = hnodes( rel(hnodes));
               neis_hrel     =   neis( rel(hnodes));
               
               krel    = length(hnodes_rel);
               knonrel = length(hnodes_nonrel);
               Wpol2 = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
               Wpol2(neis_hnonrel,:) = [zeros(knonrel,1)  Wpol(neis_hnonrel,:)] + [Wpol(hnodes_nonrel,:)  zeros(knonrel,1)];
               Wpol2(neis_hrel,:)    = [zeros(krel,1)  Wpol(neis_hrel,:)+Wpol(hnodes_rel,:)];
               Wpol = Wpol2;
           else
               Wpol(neis,:) = Wpol(neis,:) + Wpol(hnodes,:);
           end

           Es(neis) = Es(neis) - 1; %
           q = q - length(hnodes);
    
%           ind(hnodes) = 0;
           
           rel(hnodes) = [];
           n = length(rel);
           E(hnodes,:) = [];
           E(:,hnodes) = [];
           Es(hnodes)  = [];
           Wpol(hnodes,:)= [];

           hnodes = find(Es==1);
    end
end


if nnz(Es) > 0 %still have some edges
%    tmp = ECPN_C_pol(rel(ind),E(ind,ind),Wpol(ind,:));
%    tmp = ECPN_C_pol(rel,E(ind,ind),Wpol);
    tmp = ECPN_C_pol(rel,E,Wpol);
    P = poly_add(P,tmp);
end

end