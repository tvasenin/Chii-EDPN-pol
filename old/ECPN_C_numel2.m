function P = ECPN_C_numel2(VW)
%ECPN_C_numel2 Calculates ECPN for connected graph, where n = 2;
%   Detailed explanation goes here

n = length(VW);
VWpol=sym2poly_arr(VW);

P = conv(VWpol(1,:),VWpol(2,:));

syms p
P = poly2sym (P,p);
end