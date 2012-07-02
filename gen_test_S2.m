function [V E] = gen_test_S2
%gen_test_S1  Generates small test S2
%   Detailed explanation goes here

V = ones(4,1)';

E = zeros(4);
E(1,2) = 1;
E(1,3) = 1;
E(1,4) = 1;

%E(2,3) = 1;
%E(2,4) = 1;

E = (E + E');

disp('Test code:    S2');
disp('Correct ECP:  3  3  0  0');
%%strvcat('Correct answer:','[45, -17, -56, -115, -28, 408, 588, -944, -2261, 3075, 4055, -11121, 9525, -3494, 124, 264, -44, -4]')

end