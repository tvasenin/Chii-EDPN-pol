function [V E] = gen_test_chain(n)
%gen_test_chain  Generates test chain
%   Detailed explanation goes here

V = ones(n,1)';
E = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);

%ECP_correct = [1:(n-1) 0 0];


%disp(strvcat('Correct answer:', strrep(mat2str(EDP_correct),' ',', ')))

disp(['Correct ECP:  ', int2str([1:(n-1) 0 0])]);

end

