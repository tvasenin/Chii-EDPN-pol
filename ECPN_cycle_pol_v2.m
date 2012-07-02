function P = ECPN_cycle_pol_v2(rel,E,Wpol)
%ECPN_chain Calculates ECPN for a (single) cycle
%   Detailed explanation goes here
n = numel(rel);

VWpol = Wpol2VWpol(rel, Wpol);
        
switch n
    case 0
        disp('[ERROR] ECPN_cycle called with n=0')
    case 1
        disp('[ERROR] ECPN_cycle called with n=1')
    case 2
        disp('[ERROR] ECPN_cycle called with n=2')
    case 3
        %VWpol = Wpol2VWpol_opt(rel, Wpol);
        P = ECPN_full_pol(VWpol);
    case 4 %very time consuming operation on lagre graphs
        %index =  graphtraverse(E,1,'Directed',false); %traversing cycle
        %index = graphalgs('dfs',0,false,E,1,inf);% shortcut, need to place graphalgs MEX-file to the MATLAB path
        %clear E % do not need E anymore
        index = [1 0 0 0];
        index([2 4]) = find(E(1,:));
        index(3) = 10 - sum(index);
        % resorting
        rel = rel(index);
        VWpol = VWpol(index,:);
        %VWpol = Wpol2VWpol_opt(rel, Wpol(index,:));

        P = conv2(VWpol(1,:)+VWpol(3,:),VWpol(2,:)+VWpol(4,:)); %all adjacent;
        P = [0 0 P];
        
        tmp = conv2(VWpol(1,:),VWpol(3,:));
        if rel(2) || rel(4)
            P = P + [ 0 0 tmp];
        else
            P = P + [-tmp 0 0] + [0 2*tmp 0];
        end
        tmp = conv2(VWpol(2,:),VWpol(4,:));
        if rel(1) || rel(3)
            P = P + [ 0 0 tmp];
        else
            P = P + [-tmp 0 0] + [0 2*tmp 0];
        end

    otherwise
        %index =  graphtraverse(E,1,'Directed',false); %traversing cycle
        index = graphalgs('dfs',0,false,E,1,inf);% shortcut, need to place graphalgs MEX-file to the MATLAB path
        %clear E % do not need E anymore

        % resorting
        rel = rel(index);
        VWpol = VWpol(index,:);
        
        unrel_shift = ~[rel(2:end) rel(1)]';
        cum_cw = zeros(n);
        for i = 2:n
            cum_cw(:,i) = [cum_cw(2:n,i-1);  cum_cw(1,i-1)] + unrel_shift;
        end
        % rewrite for 
        
        % k = 1, i = n
        P = zeros(1,2*length(VWpol(1,:))-1);
        for i = 1:n-1 %k=1; i=1:n-1; j=i+1
            P = P + conv2(VWpol(i,:),VWpol(i+1,:));
        end
        %P = P + conv2(VWpol(n,:),VWpol(1,:)); does not count because of next loop
        for i = 1:(n-2)
            for j = (i+2):n
                k = j-i;
                tmp = conv2(VWpol(i,:),VWpol(j,:));
                cw  = cum_cw(i,k);
                ccw = cum_cw(j,n-k);
                tmp = [zeros(1,ccw) tmp zeros(1,cw)]+[zeros(1,cw) tmp zeros(1,ccw)] - [tmp zeros(1,cw+ccw)];
                P = poly_add(P,tmp);
            end
        end
        %if mod(n/2) have an odd graph
        %    c = (n/2)+1;
        
end

end