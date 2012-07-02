function P = ECPN_C_numel3(V,E,VW)
%ECPN_C_numel3 Calculates ECPN for connected graph, where n = 3;
%   Detailed explanation goes here

Es = sum(E);
if sum(Es) == 6 %full graph;
    VWpol = sym2poly_arr(VW);
    P = ECPN_full_pol(VWpol);
    syms p
    P = poly2sym(P,p);
else % 3-node chain (connected, nonfull graph  (i1)--(i2)--(i3))
    i2 = find(Es==2,1);
    i3 = mod(i2,3)+1; %to rewrite via circshift
    i1 = mod(i3,3)+1;

    P = VW(i2)*(VW(i1)+VW(i3)) + VW(i1)*VW(i3)*V(i2); %to rewrite via circshift
end

end