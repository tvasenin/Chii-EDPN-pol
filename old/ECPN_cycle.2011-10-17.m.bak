function P = ECPN_cycle(V,E,W)
%ECPN_chain Calculates ECPN for a (single) cycle
%   Detailed explanation goes here

n = numel(V);

switch n
    case 0
        disp('[ERROR] ECPN_cycle called with n=0')
    case 1
        disp('[ERROR] ECPN_cycle called with n=1')
    case 2
        disp('[ERROR] ECPN_cycle called with n=2')
    case 3
        P = ECPN_full(V.*W);
    case 4 %very time consuming operation on lagre graphs
        index =  graphtraverse(E,1,'Directed',false); %traversing cycle
        V = V(index);
        W = W(index);% resorting (we do not need E anymore)
        

% %        mult(3) = V(1)*V(2);
% %        mult(2) = V(1);
% %        mult(1) = 1;
%         %mult = V(4)*mult;
%         
%         VW = V.*W;
%         Q = 1-V;
%         %V(4) = 1; % will not use it later anyway
%         VW(4) = W(4); % CAREFUL!
%         %P = P + mult(1) * (VW(1) * W(4) + Q(1) * ECPN_ordered_chain(V(2:4),W(2:4)));
% %        P =     mult(1) * (VW(1) * W(4) + Q(1) * (VW([3 4 2])*(VW(2:4).*[1 1 V(3)])') ); %ECPN_ordered_chain(V(2:4),W(2:4)) );
% %        P = P + mult(2) * (W(4) + W(1))*(Q(2) * VW(3) + VW(2)); %prod(V(3:4).*W(3:4));
% %        P = P + mult(3) * VW(3) * (W(4) + W(1) + W(2)); % case i = n-1;
% 
%         P = V(4) *                                                          ...
%             [ (VW(1) * W(4) + Q(1) * (VW([3 4 2])*(VW(2:4).*[1 1 V(3)])') ) ... %ECPN_ordered_chain(V(2:4),W(2:4)) );
%              (W(4)+W(1))*(Q(2)*VW(3) + VW(2))                              ... %prod(V(3:4).*W(3:4));
%               VW(3) * (W([1 2 4]) * [1 1 1]')                               ] ... % case i = n-1;
%                  * [1 V(1) V(1)*V(2)]';
%         
%         %P = P * V(4);
%         
%         %P = P + Q(4) * ECPN_ordered_chain(V(1:3),W(1:3));
%         P = P + Q(4) * (VW([2 3 1])*(VW(1:3).*[1 1 V(2)])'); %it was initialization before
        P = feval(symengine, 'mu_cycle_comp1',V,W); 
        
    otherwise
        index =  graphtraverse(E,1,'Directed',false); %traversing cycle
        V = V(index);
        W = W(index);% resorting (we do not need E anymore)
        VW = V.*W;

%        P = feval(symengine, 'mu_cycle_comp2',V,W,VW);

        %need to optimize chain reductions
        %P = 0; 
%        P = P + ECPN_full(W);

        %mult = cumprod([V(n) V(1:n-1) 1]);
        mult(n+1) = sym(1); %init 
        mult(1) = 1; %init 
        for i=1:n-1
            mult(i+1) = V(i)*mult(i);
        end
        mult = V(n)*mult;
        
        Q = 1-V;
        ind_q = true(1,n-1);
        %//P = Q(n) * ECPN_ordered_chain(V(ind_q),W(ind_q),VW(ind_q)); %initialization
        P = Q(n) * feval(symengine, 'ECPN_ordered_chain_mu',V(ind_q),W(ind_q),VW(ind_q)); %initialization
        
        V(n) = 1;          
        W_n = W(n); % CAREFUL!
        for i=1:n-2
            %ind_q = i+1:n;
            ind_q = [false(1,i) true(1,n-i-1) false]; %to accomodate W_n
            %//P = P + mult(i) * Q(i) * ECPN_ordered_chain([V(ind_q) 1],[W(ind_q) W_n],[VW(ind_q) W_n]);
            P = P + mult(i) * Q(i) * feval(symengine, 'ECPN_ordered_chain_mu',[V(ind_q) 1],[W(ind_q) W_n],[VW(ind_q) W_n]);
            %V(i) = 1; not used anymore anyway
            P = P + mult(i+1) * W(i) * W_n; % contracting 2 reliable nodes
			W_n = W_n + W(i); %rewrite withiyt sym.subsref
        end
            P = P + mult(n) * W(n-1) * W_n; % case i = n-1;
        return
end

end