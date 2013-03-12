function [V E ref] = gen_test_circle(n)
%gen_test_circle  Generates test circle
%   Detailed explanation goes here

V = ones(1,n);
%E = diag(ones(n-1,1),1) + diag(ones(n-1,1),-1);
E = logical(diag(ones(n-1,1),1)) | logical(diag(ones(n-1,1),-1));

E(1,end) = 1;
E(end,1) = 1;

ref = [ -sum(2:n-2) ones(1,n-2)*n 0 0];

%disp(strvcat('Correct answer:', strrep(mat2str(EDP_correct),' ',', ')))
disp(['Correct ECP:  ', int2str(ref)]);

end

