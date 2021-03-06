function [V E] = gen_test_tree_balanced(m,h)
%gen_test_tree_balanced Generates test tree
%   Detailed explanation goes here


for i = h:-1:1
	lev_cnt(i) = m^(i-1);
end

n = sum(lev_cnt);

V = ones(1,n);

E = false(n);

lev_finish = cumsum(lev_cnt);
lev_start  = lev_finish - lev_cnt + ones(1,h);
for i = 2:h
    for j = 0:lev_cnt(h)-1
        child  = lev_start(i)   + j;
        parent = lev_start(i-1) + floor(j/m);
        E(child,parent) = 1;
    end
end


E = (E + E');


if     ([m h] == [1 1]), disp('Corr ECP 1/1 :0');
elseif ([m h] == [1 2]), disp('Corr ECP 1/2 :1  0  0');
elseif ([m h] == [1 3]), disp('Corr ECP 1/3 :');
elseif ([m h] == [2 1]), disp('Corr ECP 2/1 :0');
elseif ([m h] == [2 2]), disp('Corr ECP 2/2 :');
elseif ([m h] == [2 3]), disp('Corr ECP 2/3 :');
elseif ([m h] == [3 1]), disp('Corr ECP 3/1 :0');
elseif ([m h] == [3 2]), disp('Corr ECP 3/2 :3  3  0  0');
elseif ([m h] == [3 3]), disp('Corr ECP 3/3 :27  18  21  12   0   0');
end
    
    %EDP_correct = [sum(1:(n-1)) -(n-1):-1];

%ECP coeffs:   4  4  7  6  0  0
%disp(strvcat('Correct answer:', strrep(mat2str(EDP_correct),' ',', ')))

end

