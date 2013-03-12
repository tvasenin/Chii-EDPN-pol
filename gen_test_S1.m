function [V E ref] = gen_test_S1
%gen_test_S1  Generates small test S1
%   Detailed explanation goes here

V = ones(4,1)';

E = false(4);
E(1,2) = 1;
E(1,3) = 1;
E(1,4) = 1;

E(2,3) = 1;
E(2,4) = 1;

E = (E | E');

ref_str = '  -1  2  5  0  0';
ref = sscanf(ref_str,'%d');

disp('Test code:    S1');
disp(['Correct ECP:' ref_str]);

end