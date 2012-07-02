function Wlast = get_cum_Wlast(rel,Wpol)
%ECPN_ordered_chain Calculates ECPN for a (single) ordered chain
%   Detailed explanation goes here

n = length(rel);
switch n
    case 0
        disp('[ERROR] Attemting co call get_cum_Wlast with n=0')
%        Wlast = [];
        error('[ERROR] Attemting co call get_cum_Wlast with n=0');
    case 1
%        disp('[WARNING] Attemting co call get_cum_Wlast with n=1')
        Wlast = Wpol;
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
        
        for j = 1:n-1
            if rel(j) % [1 1]
                Wpol(j+1,:) = Wpol(j+1,:) + Wpol(j,:);
            else
                Wpol = [zeros(n,1) Wpol];
                Wpol(j+1,:) = Wpol(j+1,:) + [Wpol(j,2:end) 0];
            end
        end
        Wlast = Wpol(end,:);
end
end