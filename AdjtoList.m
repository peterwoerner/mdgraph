function [l] = AdjtoList(A)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
[n,~] = size(A);
l = [];
for i = 1:n
    for j = i:n
        if A(i,j) ~= 0
            l = [l; i, j, A(i,j)];
        end
    end
end
            

end

