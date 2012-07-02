function Q = rel2Qpol(rel)
%Q becomes logical

Q(:,1) = -~rel;
Q(:,2) = -rel + 1;
 
end