function B = sym2poly_arr(A)

n = numel(A);

B=0;
for i=n:-1:1
    tmp = fliplr(sym2poly(A(i)));
    B(i,1:length(tmp)) = tmp;
end
B = fliplr(B);

end