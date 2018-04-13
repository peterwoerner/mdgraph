function [rij, dir] = distancematrix(r1, r2, varargin)
%distance matrix Creates a matrix of the distances between points in vectors r1, and r2 and unit direction vectors
%   Detailed explanation goes here
[na, dim] = size(r1);
[na2, ~] = size(r2);
for i = 1:na
    for j = 1:na2
        dr = 0;
        for k = 1:dim
            dr = dr + (r1(i,k)-r2(j,k))^2;
        end
        dr = sqrt(dr);
        rij(i,j) = dr;
        for k = 1:dim
            if dr ~= 0
                dir(i,j,k) = (r1(i,k)-r2(j,k))/dr;
            else
                dir(i,j,k) = 0;
            end
        end
    end
end
end

