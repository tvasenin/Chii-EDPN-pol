function V = rel2V(rel)
%V becomes logical
syms p
V = rel + p * (~rel);
%V(1:length(rel)) = p;
%V(rel) = 1;
end