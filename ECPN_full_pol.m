function P = ECPN_full_pol(VWpol)
%ECPN_full Calculates ECPN for full G(VW)
%   Detailed explanation goes here

n = length(VWpol(:,1));
switch n
    case 0
        disp('[ERROR] Attempting to call ECPN_full with n=0')
    case 1
        disp('[ERROR] Attempting to call ECPN_full with n=1')
    case 2
        %disp('[WARNING] Attempt to call ECPN_full with n=2')
%        P = ECPN_C_numel2_pol(VW);
        P = conv2(VWpol(1,:),VWpol(2,:));
    case 3
%        disp('[WARNING] Attempt to call ECPN_full with n=3')
%        P = sum(VW.*circshift(VW,[0 1]));
%        P = VW([2 3 1])*VW';
        P = conv2(VWpol(1,:)+VWpol(2,:),VWpol(3,:)) + conv2(VWpol(1,:),VWpol(2,:));
    case 4
        P = conv2(VWpol(1,:)+VWpol(2,:), VWpol(3,:)+VWpol(4,:)) + conv2(VWpol(1,:), VWpol(2,:)) + conv2(VWpol(3,:),VWpol(4,:));
    case 5 
        P = conv2(VWpol(1,:)+VWpol(2,:), VWpol(3,:)+VWpol(4,:)+VWpol(5,:)) + conv2(VWpol(3,:)+VWpol(4,:), VWpol(5,:)) + conv2(VWpol(1,:),VWpol(2,:)) + conv2(VWpol(3,:),VWpol(4,:));
    otherwise
  
%        P = sum(sum(triu(VW'*VW,1))); %profiling shows it's the best we can do
%        P = ones(1,n) * triu(VW'*VW,1) * ones(n,1); %well, it's FAR better :)
        
        %v1
        P(2*numel(VWpol(1,:))-1) = 0;
        for i = 1:(n-1)
            P = P + conv2(VWpol(i,:),sum(VWpol(i+1:n,:),1));
        end
        
        %v2 (wrong at e.g. 300 nodes)
        %i1 = ceil(n/2);
        %i2 = floor(n/2);
        %P = conv2(sum(VWpol(i1:n,:),1),sum(VWpol(1:i2,:),1));
        %P = P + ECPN_full_pol(VWpol(1:i1-1,:)) + ECPN_full_pol(VWpol(i2+1:n,:));
end

end