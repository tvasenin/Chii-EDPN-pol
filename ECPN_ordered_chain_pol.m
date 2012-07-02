function P = ECPN_ordered_chain_pol(rel,VW)
%ECPN_ordered_chain Calculates ECPN for a (single) ordered chain
%   Detailed explanation goes here

n = length(rel);
switch n
    case 0
        disp('[WARNING] Attemting co call ECPN_ordered_chain with n=0')
%        P = 0;
    case 1
%        disp('[WARNING] Attemting co call ECPN_ordered_chain with n=1')
        P = 0;
    case 2
        %P = ECPN_C_numel2_pol(VW);
        P = conv2(VW(1,:),VW(2,:)); % need for speed
    case 3
%        VW = V.*W;
%        P = ECPN_C_numel3 (V,E,V.*W);
%        P = VW(2)*(VW(1)+VW(3)) + VW(1)*VW(3)*V(2);
%        P = VW([2 3 1])*(VW.*[1 1 V(2)])';
        if rel(2) 
            P = conv2(VW(1,:)+VW(3,:),VW(2,:)) + conv2(VW(3,:),VW(1,:));
        else
            P = [0 conv2(VW(1,:)+VW(3,:),VW(2,:))] + [conv2(VW(3,:),VW(1,:)) 0];
        end
    otherwise
        %for i = 2:n-1 %may need optimization
        %    W(i) = W(i) + V(i-1)*W(i-1); %can't use VW(i-1)
        %end
        
        %VVW = V.*circshift(VW,[0,-1]); %VVW(i) = V(i) * VW(i+1);
        %VVW(n) = 0;
        %P = VVW * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));

        %VVW = V.*[VW(2:n) 0]; %VVW(i) = V(i) * VW(i+1);
        %P = VVW * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));
        %P = (V.*[VW(2:n) 0]) * W'; %P = sum(VVW(1:n-1) .* W(1:n-1));
        
        VW1 = VW(1,:);
        P = conv2(VW(2,:),VW1);
        %cumnrel = cumsum(~rel);
        for i = 3:n
            if ~rel(i-1)
            %    VW  = padarray(VW, [0,1],'pre');
                VW1 = [VW1 0]; %padarray(VWn,[0,1],'post');
            %    P   = [0 0 P]; %padarray(P,  [0,2],'pre');
            end
            %VWn = VW(i-1,:) + VWn;
            VW1 = poly_add(VW1,VW(i-1,:));
            %P = P + conv2(VW(i,:),VW1);
            P = poly_add(P, conv2(VW(i,:),VW1));
        end
end
%P = polytrim_fast(P);
end