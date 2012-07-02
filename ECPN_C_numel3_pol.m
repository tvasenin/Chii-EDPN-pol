function P = ECPN_C_numel3_pol(rel,E,VW)
%ECPN_C_numel3 Calculates ECPN for connected graph, where n = 3;
%   Detailed explanation goes here

Es = sum(E);
if sum(Es) == 6 %full graph;
%    P = ECPN_full_pol(VW);
    P = conv2(VW(1,:)+VW(2,:),VW(3,:)) + conv2(VW(1,:),VW(2,:));
else % 3-node chain (connected, nonfull graph  (i1)--(i2)--(i3))
    i2 = find(Es==2,1);
    i3 = mod(i2,3)+1; %to rewrite via circshift
    i1 = mod(i3,3)+1;

    if rel(i2) % V(i2) == 1
        P = conv2(VW(i1,:)+VW(i3,:),VW(i2,:)) + conv2(VW(i3,:),VW(i1,:));
    else       % V(i2) == p
        %P = VW(i2)*(VW(i1,:)+VW(i3,:)) + VW(i1)*VW(i3)*V(i2); %to rewrite via circshift
        P = [0 conv2(VW(i1,:)+VW(i3,:),VW(i2,:))] + [conv2(VW(i3,:),VW(i1,:)) 0];
    end
end

end