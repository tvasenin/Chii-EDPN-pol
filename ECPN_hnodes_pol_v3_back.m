function [P, rel, E, Wpol] = ECPN_hnodes_pol_v3_back(rel,E,Wpol)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here

n = numel(rel);
Es = sum(E);
P = zeros(1,2*length(Wpol(1,:))+1);

hnodes = find(Es==1); %maybe it's possible to rewrite via hnodes = Es==1;
%ind = true(1,n);
%q = sum(Es)/2;
%disp('===================================')
while ~isempty(hnodes)
    switch n  %q
        case 2  % have only 2 nodes
            %P = P + ECPN_C_numel2(VW(hnodes)); % now we recalc VW.
            P = poly_add(P,ECPN_C_numel2_pol_v2(rel,Wpol)); % now we recalc VW.
            E = [];
            rel = [];
            Wpol = [];
            break
        otherwise

            if length(hnodes) > 1 % actally need standalone code path here
                [i1, neis] = find(E(hnodes,:));
                [~, i2] = sort(i1); %to optimize!
                neis = neis(i2);
                
                [neis, ind2] = unique_fast(neis); % to optimize
                hnodes = hnodes(ind2);
            else % have only one hnode
                neis = find(E(hnodes,:),1);
            end
    
%           if Es(hnei) == 2
%               disp('WOW! Hanging chain has been found!');
%           end
           %P = P + VW(hnodes)*VW(hneis)'; %shouldn't use VW for hnode!!
           %tmp = zeros(1,length(hnodes));
           
           for i = 1:length(hnodes)
               hnode = hnodes(i);
               nei = neis(i);
               if      rel(hnode) &&  rel(nei) % [1 1]
                   P = P + [0  0  conv2(Wpol(hnode,:),Wpol(nei,:))];
%                   P(3:end) = P(3:end) + conv2(Wpol(hnode,:),Wpol(nei,:));
               elseif ~rel(hnode) && ~rel(nei) % [p p]
                   P = P + [conv2(Wpol(hnode,:),Wpol(nei,:))  0  0];
%                   P(1:end-2) = P(1:end-2) + conv2(Wpol(hnode,:),Wpol(nei,:));
               else                            % [0 p] or [p 0]
                   P = P + [0  conv2(Wpol(hnode,:),Wpol(nei,:))  0];
%                   P(2:end-1) = P(2:end-1) + conv2(Wpol(hnode,:),Wpol(nei,:));
               end
           end
           
           %W(hneis) = W(hneis) + VW(hnodes);
           %Wpol = padarray (Wpol,[0,1],'pre'); % now dim(Wpol) should be equal to dim(VWpol)
           
           %Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
           %Wpol(neis,:) = Wpol(neis,:) + VWpol(hnodes,:);
               
           if length(hnodes) > 1
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
                   P = [0 0 P];
               else
                   Wpol(neis,:) = Wpol(neis,:) + Wpol(hnodes,:);
               end
           else %have only one hnode, can use hnode and nei
               if rel(hnode)
                   Wpol(nei,:) = Wpol(nei,:) + Wpol(hnode,:);
               else
                   Wpol2 = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
                   Wpol2(nei,:) = [ 0  Wpol(nei,:) ] + [ Wpol(hnode,:)  0 ];
                   Wpol = Wpol2;
                   P = [0 0 P];
               end
           end

           Es(neis) = Es(neis) - 1; %
%           q = q - length(hnodes);
    
           rel(hnodes) = [];
           n = length(rel);
           E(hnodes,:) = [];
           E(:,hnodes) = [];
           Es(hnodes)  = [];
           Wpol(hnodes,:)= [];

           hnodes = find(Es==1);
    end
end

%P = P(find(P,1):end); % trim leading zeros
%P = polytrim_fast(P); % trim leading zeros

end