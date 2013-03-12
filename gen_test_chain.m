function [V E ref] = gen_test_chain(n)
%gen_test_chain  Generates test chain
%   Detailed explanation goes here

V = ones(1,n);
%E = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
E = logical(diag(ones(n-1,1),1)) | logical(diag(ones(n-1,1),-1));

ref = [1:(n-1) 0 0];

disp(['Correct ECP:  ', int2str(ref)]);

end

