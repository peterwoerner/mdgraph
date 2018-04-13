function [Ar, Af, Au] = getA(r, q, epsilon)
%Gets adjacency matrices of network of force and network of potential energy
%assuming a Lennard Jones potential.  r is the positions of the atoms,
%epsilon is the depth of the potential well and sigma is the finite
%distance at which the inter-particle potential is zero.
[l,n] = size(r);  %l is the length of r (number of atoms), n is the dimension (i.e number of degrees of freedom in the position)
Af = zeros(l,l); % force matrix
Au = zeros(l,l); %potential energy matrix
Ar = zeros(l,l); %distance matrix
[Ar,~] = distancematrix(r, r); 
q = reshape(q, [l,1]);
Au = (q*transpose(q))./(Ar)/epsilon;
Af = (q*transpose(q))./(Ar.*Ar)/epsilon;
for i = 1:l
   Au(i,i) = 0;
   Af(i,i) = 0;
   Ar(i,i) = 0;
end

end

