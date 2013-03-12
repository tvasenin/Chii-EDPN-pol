function [ c ] = poly_add_safe( a, b )
%poly_add Adds 2 polynoms
%   Detailed explanation goes here

k = length(b) - length(a);

%c = [zeros(1,k) a] + [zeros(1,-k) b];
if (k > 0)
c = [zeros(1,k)  a] + b;
%c = [b(1, k) (b( (k+1):end) + a)];
else
c = [zeros(1,-k) b] + a;
%c = [a(1,-k) (a((-k+1):end) + b)];
end

end
