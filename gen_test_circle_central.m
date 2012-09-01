function [V E ref] = gen_test_circle_central(n)
%gen_test_circle  Generates test circle
%   Detailed explanation goes here

[V E ref] = gen_test_circle(n);

central_perm = [];
for i = 1:n
    if mod(i,2)
       central_perm = [central_perm i];
    else
       central_perm = [i central_perm];
    end
end

central_perm2 = zeros(1,n);

for i = 1:n
    central_perm2(i) = find(central_perm == i);
end


E = E(central_perm2,central_perm2);

%EDP_correct = [sum(1:(n-1)) -(n-1):-1];

%disp(strvcat('Correct answer:', strrep(mat2str(EDP_correct),' ',', ')))

end

