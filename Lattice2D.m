function [r] = Lattice2D(a, siz, type)
%

periodicity = true;
sq2 = sqrt(2);
if strcmp(type, 'sq2') || strcmp(type, 'hex')    
  ri = [0,0; 1/2, 1/2]; 
elseif strcmp(type, 'sq')
  ri = [0,0];
else
  error('Type selected invalid. Only types "sq" and "sq2" are available')
end  

r0 = ri;
[n1, n2] = size(ri);
for i = 0:siz(1)-1
    for j = 0:siz(2)-1

        radd = ri + [ones(n1,1)*i ones(n1,1)*j];
        r0 = [r0;radd];
    end
end
r = unique(r0, 'rows');
if periodicity ==true
    
    epsilon = .1;
    %remove top, right and front boundaries (due to periodicity of the crystal)
    list = [];
    for i = 1:length(r)
        if r(i,1) > siz(1)-epsilon
            list = [list,i];
        elseif  r(i,2) > siz(2)-epsilon   
            list = [list,i];
        end
    end
    r(list,:) = [];
end

r = r*a;



end

