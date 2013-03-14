function [P, rel, E, Wpol] = ECPN_hnodes_pol_v4_back(rel,E,Wpol)
%ECPN_hnodes Calculates ECPN for a graph with hanging trees (or just chains or nodes)
%   Detailed explanation goes here

n = numel(rel);
P = zeros(1,2*length(Wpol(1,:))+1);

[lmat pmat mcnt] = GetHnodesQueue(E);

leftover = true(1,n);

for k = 1:length(mcnt)
    lhnodes = mcnt(k);
    hnodes = lmat(k,1:lhnodes);
    neis   = pmat(k,1:lhnodes);
    %P = P + VW(hnodes)*VW(hneis)'; %shouldn't use VW for hnode!!

    for i = 1:lhnodes
        hnode = hnodes(i);
        nei = neis(i);
        nrh = ~rel(hnode);
        nrn = ~rel(nei);
        tmp = conv2(Wpol(hnode,:),Wpol(nei,:));

        if     nrh && nrn  % [p p]
            P = P + [ tmp 0  0];
        elseif        nrn  % [1 p]
            P = P + [0  tmp  0];
        elseif nrh         % [p 1]
            P = P + [0  tmp  0];
        else               % [1 1]
            P = P + [0  0  tmp];
        end
    end

    %W(hneis) = W(hneis) + VW(hnodes);
    %Wpol = padarray (Wpol,[0,1],'pre'); % now dim(Wpol) should be equal to dim(VWpol)
    %Wpol = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
    %Wpol(neis,:) = Wpol(neis,:) + VWpol(hnodes,:);

    if (lhnodes == 1) %have only one hnode, can use hnode and nei
        if rel(hnode)
            Wpol(nei,:) = Wpol(nei,:) + Wpol(hnode,:);
        else
            if Wpol(hnode,1)
                Wpol2 = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
                Wpol2(nei,:) = [ 0  Wpol(nei,:) ] + [ Wpol(hnode,:)  0 ];
                Wpol = Wpol2;
                P = [0 0 P];
            else
                Wpol(nei,:) = Wpol(nei,:) + [ Wpol(hnode,2:end)  0 ];
            end
        end
    else
        hnodes_nonrel = hnodes(~rel(hnodes));
        if ~isempty(hnodes_nonrel)
            neis_hnonrel  =   neis(~rel(hnodes));
            hnodes_rel    = hnodes( rel(hnodes));
            neis_hrel     =   neis( rel(hnodes));

            krel    = length(hnodes_rel);
            knonrel = length(hnodes_nonrel);
            if find(Wpol(hnodes_nonrel,1))
                Wpol2 = [zeros(n,1) Wpol]; % padarray is sooo sloow ^^
                Wpol2(neis_hnonrel,:) = [zeros(knonrel,1)  Wpol(neis_hnonrel,:)] + [Wpol(hnodes_nonrel,:)  zeros(knonrel,1)];
                Wpol2(neis_hrel,:)    = [zeros(krel,1)  Wpol(neis_hrel,:)+Wpol(hnodes_rel,:)];
                Wpol = Wpol2;
                P = [0 0 P];
            else
                Wpol(neis_hnonrel,:) = Wpol(neis_hnonrel,:) + [Wpol(hnodes_nonrel,2:end)  zeros(knonrel,1)];
                Wpol(neis_hrel,:)    = Wpol(neis_hrel,:)+Wpol(hnodes_rel,:);
            end
        else
            Wpol(neis,:) = Wpol(neis,:) + Wpol(hnodes,:);
        end
    end
    leftover(hnodes) = false;
end

%[~, fnzc] = find(Wpol,1); % approximate
rel  =  rel(leftover);
E    =    E(leftover,leftover);
%Wpol = Wpol(leftover,fnzc:end);
Wpol = Wpol(leftover,:);


%P = P(find(P,1):end); % trim leading zeros
%P = polytrim_fast(P); % trim leading zeros

end