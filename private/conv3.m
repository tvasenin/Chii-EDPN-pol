function res = conv3(a, b)

if (a(1) == 0)
    a = polytrim_fast(a);
end
if (b(1) == 0)
    b = polytrim_fast(b);
end

res = conv2(a(find(a,1):end),b(find(b,1):end));
%res = conv2(a,b);

end
