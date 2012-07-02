function E = ContractNodes(E,nodes)
%ContractNode Contracts edges matrix by node k without deleting it
%   Detailed explanation goes here

%need optimizing
for k = nodes
    cont = E(k,:)~=0;
    E(cont,cont) = 1;
end
for i = 1:length(E)
    E(i,i) = 0;
end
%    E = spdiags(sparse(length(E),1),0,E); %E = A - diag(diag(A));
end

