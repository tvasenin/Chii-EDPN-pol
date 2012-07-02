function [V E] = gen_test_circle(n)
%gen_test_chain  Generates test circle
%   Detailed explanation goes here

V = ones(n,1)';
E = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);

E(1,end) = 1;
E(end,1) = 1;

ECP_correct = [ -sum(2:n-2) ones(1,n-2)*n 0 0];



%disp(strvcat('Correct answer:', strrep(mat2str(EDP_correct),' ',', ')))
disp(['Correct ECP:  ', int2str(ECP_correct)]);

end

