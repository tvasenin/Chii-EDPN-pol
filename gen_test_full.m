function [V E ref] = gen_test_full(n)
%gen_test_full  Generates full graph test
%   Detailed explanation goes here

V = ones(n,1)';
E = ones(n) - eye(n);

ref = [n*(n-1)/2 0 0];

disp(['Test code:    F' int2str(n)]);
disp(['Correct ECP:  ' int2str([n*(n-1)/2 0 0])]);

end
