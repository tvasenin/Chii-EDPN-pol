function P = ECPN_ordered_chain(V,W,VW)
%ECPN_ordered_chain Calculates ECPN for a (single) ordered chain
%   Detailed explanation goes here

n = numel(V);

switch n
    case 0
        disp('[WARNING] Attemting co call ECPN_ordered_chain with n=0')
%        P = 0;
    case 1
        disp('[WARNING] Attemting co call ECPN_ordered_chain with n=1')
%        P = 0;
    case 2
%        P = ECPN_C_numel2(V.*W);
        P = ECPN_C_numel2(VW);
    case 3
%        VW = V.*W;
        P = VW([2 3 1])*(VW.*[1 1 V(2)])';
%        P = VW(2)*(VW(1)+VW(3)) + VW(1)*VW(3)*V(2);
        %P = ECPN_C_numel3(V,E,V.*W);
    otherwise
       
        for i = 2:n-1 %may need optimization
            W(i) = W(i) + V(i-1)*W(i-1); %can't use VW(i-1)
        end
        
        %VVW = V.*circshift(VW,[0,-1]); %VVW(i) = V(i) * VW(i+1);
        %VVW(n) = 0;
        %P = VVW * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));

        %VVW = V.*[VW(2:n) 0]; %VVW(i) = V(i) * VW(i+1);
        %P = VVW * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));
        P = (V.*[VW(2:n) 0]) * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));

end