function [V E ref] = gen_test_big_G1
%gen_test_full  Generates big test G1
%   Detailed explanation goes here

V = ones(1,10);

E = false(10);
E(1,3) = 1;

E(2,3) = 1;
E(2,5) = 1;
E(2,7) = 1;
E(2,8) = 1;

E(3,4) = 1;
E(3,7) = 1;
E(3,8) = 1;
E(3,9) = 1;

E(4,6) = 1;
E(4,8) = 1;
E(4,9) = 1;

E(5,7) = 1;

E(6,9) = 1;

E(7,8) = 1;

E(8,9) = 1;
E(8,10)= 1;

E = (E | E');

ref_str = '  -1   6  -8 -12  17  26  17   0   0';
ref = sscanf(ref_str,'%d');

disp('Test code:    G1');
disp(['Correct ECP:' ref_str]);

end