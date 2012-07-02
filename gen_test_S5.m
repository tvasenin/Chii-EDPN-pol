function [V E] = gen_test_S5
%gen_test_S1  Generates small test S5
%   Detailed explanation goes here

syms p;

V = ones(6,1)';

E = zeros(6);
E(1,2) = 1;
E(1,4) = 1;
E(2,3) = 1;
E(2,6) = 1;
E(3,4) = 1;
E(4,5) = 1;
E(4,6) = 1;


E = (E + E');

disp('Test code:    S5');
disp('Correct ECP:  1  -2  -3  12   7   0   0');

end