function [ G ] = KaoFo2std( Kao, Fo )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

n = length(Kao)-1;
G = sparse(zeros(n));

for i = 1:n
    G(i,Fo((Kao(i)+1):(Kao(i+1)))) = 1;
end
end

