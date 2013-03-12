function [V E] = gen_test_big_G1_easy
%gen_test_full  Generates big test G1_easy
%   Detailed explanation goes here

V = ones(10,1)';

E = false(10);
E(1,3) = 1;

E(2,3) = 1;
%E(2,5) = 1;
%E(2,7) = 1;
E(2,8) = 1;

E(3,4) = 1;
E(3,7) = 1;
E(3,8) = 1;
E(3,9) = 1;

%E(4,6) = 1;
E(4,8) = 1;
E(4,9) = 1;

E(5,7) = 1;

%E(6,9) = 1;

E(7,8) = 1;

E(8,9) = 1;
E(8,10)= 1;

E = (E | E');

disp(strvcat('Correct answer with 4 commented edges:','[45, -13, -40, -63, 40, 232, -148, -258, 239, 53, -95, 11, 6]'))

end