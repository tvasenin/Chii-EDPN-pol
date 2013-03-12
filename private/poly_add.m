function [ a ] = poly_add( a, b )
%poly_add Adds 2 polynoms
%   Detailed explanation goes here

%if (a(1) == 0)
%    a = polytrim_fast(a);
%end
%if (b(1) == 0)
%    b = polytrim_fast(b);
%end

k = length(b) - length(a);
a = [zeros(1,k) a] + [zeros(1,-k) b];

end

