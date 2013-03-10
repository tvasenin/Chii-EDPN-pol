function [ a ] = poly_add( a, b )
%poly_add Adds 2 polynoms
%   Detailed explanation goes here

k = length(b) - length(a);
a = [zeros(1,k) a] + [zeros(1,-k) b];

end

