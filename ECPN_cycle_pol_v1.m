function P = ECPN_cycle_pol_v1(rel,E,Wpol)
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
        %VWpol = Wpol2VWpol(rel, Wpol);
        %index =  graphtraverse(E,1,'Directed',false); %traversing cycle
        index = graphalgs('dfs',0,false,E,1,inf);% shortcut, need to place graphalgs MEX-file to the MATLAB path
        %P = feval(symengine, 'muCycleComp',V(index),W(index));%,(V.*W));
        %clear E % do not need E anymore

        % resorting
        rel = rel(index);
        VWpol = VWpol(index,:);
        Wpol = Wpol(index,:);
        
        
        %need to optimize chain reductions
% %        P = P + ECPN_full(W);
% 
        %Qpol = rel2Qpol(rel); % do not need it for manual convolution
        %mult = cumprod([V(n) V(1:n-1) 1]);
        %mult2 = V(n);
        cum_nonrel = cumsum(~rel([n,1:(n-1)]));
        
        P = 0;
        if ~rel(n)
            rel(n) = 1; % outside of first ind_q
            ind_q = true(1,n-1);
            %P = Q(n) * feval(symengine, 'ECPN_ordered_chain_mu',V(ind_q),VW(ind_q)); %initialization
            %P = conv( Qpol(n,:), ECPN_ordered_chain_pol(rel(ind_q),VWpol(ind_q,:)) );
            %manual convolution
            P = ECPN_ordered_chain_pol(rel(ind_q),VWpol(ind_q,:));
            P = [-P 0] + [0 P];
        end
            
        %VW(n) = W(n); % CAREFUL! but it's outside of ind_q
        VWpol(n,:) = [0 Wpol(n,:)];
    

        for i=1:n-2
            if ~rel(i)
%             ind_q = [false(1,i) true(1,n-i-1) false]; %to accomodate W_n
                ind_q = i+1:n;
%                P = P + mult(i) * Q(i) * ECPN_ordered_chain([V(ind_q) 1],[VW(ind_q) W_n]);
%                P = P + mult2 * Q(i) * feval(symengine, 'ECPN_ordered_chain_mu',[V(ind_q)],[VW(ind_q)]);
                tmp = ECPN_ordered_chain_pol(rel(ind_q),VWpol(ind_q,:));
                %tmp = conv(Qpol(i,:),tmp);
                tmp = [-tmp 0] + [0 tmp];
                P = poly_add(P,[tmp  zeros(1,cum_nonrel(i))]);
            end
        
%             %V(i) = 1; not used anymore anyway
%             P = P + mult(i+1) * W(i) * W_n; % contracting 2 reliable nodes
            %P = P + mult2 * VW(i) * VW(n); % contracting 2 reliable nodes
            
            tmp = [conv2(VWpol(i,:),VWpol(n,:)) zeros(1,cum_nonrel(i))];
            P = poly_add(P,tmp);
%           %P =  P + tmp;
% 			W_n = W_n + W(i); %rewrite withiyt sym.subsref
            %VW(n) = VW(n) + W(i);
            %tmp = poly_add(VWpol(n,:),Wpol(i,:));
            %VWpol(n,1:length(tmp)) = tmp;
            VWpol(n,:) = VWpol(n,:) + [0 Wpol(i,:)];
            %mult2 = mult2 * V(i);
        
        end
        %P = P + mult2 * VW(n-1) * VW(n); % case i = n-1;
        tmp = [conv2(VWpol(n-1,:),VWpol(n,:)) zeros(1,cum_nonrel(n-1))];
        P = poly_add(P,tmp);
end

end