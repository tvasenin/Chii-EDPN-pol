function P = ECPN_C_nCH_nCC_numel4(V,E,W);
%ECPN_C_nCH_nCC_numel4 Calculates ECPN for connected graph, where n = 4;
%   Detailed explanation goes here

n = numel(V);
Es = sum(E);
q  = sum(Es)/2;

if q == n*(n-1)/2
    P = ECPN_full(V.*W);
    return
end 

[Es, index] = sort(Es,'descend'); % sorting nodes by degree
V = V(index);
W = W(index);


VW = V.*W;

P = V(1)*( W(1)*sum(VW(2:4)) + sum(VW(2:4).*VW([3 4 2])) ); %circshift
if q == 3  %Es(2) == 1 %Es == [3 1 1 1]
elseif q == 4  %Es(2) == 2 %Es == [3 2 2 1]
    %P = prod(VW([2 3]))  +  V(1)*( W(1)*sum(VW([2 3 4])) + sum(VW([3 4]).*VW([4 2])) );
    P = P + (1-V(1)) * VW(2)*VW(3);
end
if q == 5 %Es(2) == 3 %Es == [3 3 2 2]
    P = P + (1-V(1)) * ( VW(2)*(VW(3)+VW(4)) + V(2)*VW(3)*VW(4) ); %to rewrite via circshift
end


end