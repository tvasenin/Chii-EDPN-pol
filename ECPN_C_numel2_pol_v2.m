function P = ECPN_C_numel2_pol_v2(rel,Wpol)
%ECPN_C_numel2 Calculates ECPN for connected graph, where n = 2;
%   Detailed explanation goes here

if rel(1)
    if rel(2)
        P = conv2(Wpol(1,:),Wpol(2,:));
    else
        P = conv2(Wpol(1,:),[Wpol(2,:) 0]);
    end
elseif rel(2)
    P = conv2([Wpol(1,:) 0],Wpol(2,:));
else
    P = [conv2(Wpol(1,:),Wpol(2,:)) 0];
end

end