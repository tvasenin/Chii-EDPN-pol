function P = ECPN_cycle(V,E,W)
%ECPN_chain Calculates ECPN for a (single) cycle
%   Detailed explanation goes here

n = numel(V);

switch n
%    case 1
%    case 2
    case 3
        P = ECPN_C_numel3(V,E,V.*W);
    otherwise % including n = 4
        index =  graphtraverse(E,1,'Directed',false); %traversing cycle from first node
        V = V(index); W = W(index);% resorting (we do not need E anymore)

        %dynamic precalc of symbolic cumprod array
        cp = [V ; zeros(n-1,n)];
        for i = 1:n-1
            cp(i+1,1:n-i) = cp(i,1:n-i).*V(i+1:n);
        end
       
        %cp = [V ; zeros(n-1,n)];
        %for i = 1:n-1
        %    cp(i+1,:) = cp(i,:).*circshift(V,[0,-i]);
        %end
        %cp = cp';
        
        %need to heavily optimizie!
        WW = W'*W;
        P = 0;

%        for i = 1:n-1
%            P = P + WW(i,i+1:n) .* (cp(i,2:n-i+1) + cp(j,n-j+i+1));
%            % for j = i+1:n
%            %     P = P + WW(i,j) * (cp(i,j-i+1) + cp(j,n-j+i+1));
%            % end
%        end
%        P = P - cp(1,n)*sum(sum(triu(WW,1)));
       
        for i = n:-1:1
            cp_antidiag(i) = cp(n+1-i,i);
        end
        cp_antidiag = cp_antidiag';
        
        for i = 1:n-1
%            P = P + sum(WW(i+1:n,i).*(cp(2:n-i+1,i) + cp(i,1)*cp_antidiag(i+1:n)));
            ololo = WW(i+1:n,i).*(cp(2:n-i+1,i) + cp(i,1)*cp_antidiag(i+1:n));
            P = P + sum(ololo);
        end
        P = P - cp(n,1)*sum(sum(triu(WW,1)));
end

end