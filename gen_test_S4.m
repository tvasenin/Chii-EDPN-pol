function [V E] = gen_test_S4
%gen_test_S1  Generates small test S4
%   Detailed explanation goes here

syms p;

V = ones(5,1)';


E = zeros(5);
E(1,3) = 1;
E(1,4) = 1;
E(1,5) = 1;
E(2,3) = 1;
E(2,4) = 1;
E(2,5) = 1;

E = (E + E');

disp('Test code:    S4');
disp('Correct ECP:  1 -6  9  6  0  0');

end