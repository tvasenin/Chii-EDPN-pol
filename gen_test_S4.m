function [V E ref] = gen_test_S4
%gen_test_S4  Generates small test S4
%   Detailed explanation goes here

V = ones(5,1)';

E = false(5);
E(1,3) = 1;
E(1,4) = 1;
E(1,5) = 1;
E(2,3) = 1;
E(2,4) = 1;
E(2,5) = 1;

E = (E | E');

ref_str = '  1 -6  9  6  0  0';
ref = sscanf(ref_str,'%d');

disp('Test code:    S4');
disp(['Correct ECP:' ref_str]);

end