function [ r ] = GenerateHole(r, location, radius )
%#Pass in atom positions hole location and hole radius.  All atoms with in hole will be deleted.
delete_list = [];
[n1, n2] = size(r);
for i = 1:n1
    if sqrt((r(i,1)-location(1))^2 + (r(i,2)-location(2))^2) < radius
        delete_list = [delete_list, i];
    end
end
r(delete_list,:) = [];



end

