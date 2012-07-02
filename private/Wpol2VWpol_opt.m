function [ VWpol ] = Wpol2VWpol_opt( rel, Wpol )
%Wpol2VWpol_opt Apply rel to Wpol but without excessive zero padding
%   Detailed explanation goes here

if nnz(rel)
    for i = numel(rel):-1:1
        if rel(i)
            VWpol(i,:) = [0 Wpol(i,:)];
        else
            VWpol(i,:) = [Wpol(i,:),0];
        end
    end
else
    VWpol = Wpol;
end
end
