function P = ECPN_C_numel2_pol(VW)
%ECPN_C_numel2 Calculates ECPN for connected graph, where n = 2;
%   Detailed explanation goes here

P = conv2(VW(1,:),VW(2,:));

end