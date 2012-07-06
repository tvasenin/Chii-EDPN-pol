function P = ECPN_C_numel2_pol_v2(rel,Wpol)
%ECPN_C_numel2 Calculates ECPN for connected graph, where n = 2;
%   Detailed explanation goes here

% P has various resulting size in these four cases
if rel(1)
    if rel(2) % [1 1]
        P = conv2(Wpol(1,:),Wpol(2,:));
    else      % [1 0]
        P = conv2(Wpol(1,:),[Wpol(2,:) 0]);
    end
elseif rel(2) % [0 1]
    P = conv2([Wpol(1,:) 0],Wpol(2,:));
else          % [0 0]
    P = [conv2(Wpol(1,:),Wpol(2,:)) 0 0];
end

end