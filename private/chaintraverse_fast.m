function index = chaintraverse_fast(E)
%CHAINTRAVERSE Traverses the chain from starting node
%   Detailed explanation goes here

n = length(E); %assert n>2 !!!
%index = zeros(1,n);
%index([1 n]) = find(sum(E)==1);
%i1 = index(1);
%i2 = find(E(i1,:));
%index(2) = i2;

%index = zeros(1,n);
i3 = find(sum(E)==1);
i1 = i3(1);
i2 = i3(2);
index(n) = i2;
index(1) = i1;
i2 = find(E(i1,:));
index(2) = i2;


%[c, r] = find(tril(E)); % can use [c r] instead of [r c] because E is symmetric

for i = 3:n-1
    i3 = find(E(i2,:));
    if i3(1) == i1
        i1 = i2;
        i2 = i3(2);
    else
        i1 = i2;
        i2 = i3(1);
    end
    index(i) = i2;
end

end
