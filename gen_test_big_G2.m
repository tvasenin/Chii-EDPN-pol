function [V E ref] = gen_test_big_G2
%gen_test_full  Generates big test G2
%   Detailed explanation goes here

V = ones(1,10);

E = false(10);
E(1,2) = 1;
E(1,4) = 1;

E(2,3) = 1;
E(2,4) = 1;

E(3,6) = 1;
E(3,7) = 1;

E(4,5) = 1;

E(5,6) = 1;
E(5,8) = 1;
E(5,9) = 1;
E(5,10)= 1;

E(6,7) = 1;
E(6,9) = 1;
E(6,10)= 1;

E(7,10)= 1;
E(8,9) = 1;
E(9,10)= 1;

E = (E | E');

ref_str = '  4   5 -13 -22   9  23  22  17   0   0';
ref = sscanf(ref_str,'%d');

disp('Test code:    G2');
disp(['Correct ECP:' ref_str]);

end