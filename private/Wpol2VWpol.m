function [ VWpol ] = Wpol2VWpol( rel, Wpol )
%Wpol2VWpol Apply rel to Wpol
%   Detailed explanation goes here

%VWpol = zeros(numel(rel),numel(Wpol(1,:))+1);
for i = numel(rel):-1:1
    if rel(i)
        VWpol(i,:) = [0 Wpol(i,:)];
    else
        VWpol(i,:) = [Wpol(i,:),0];
    end
end

end

