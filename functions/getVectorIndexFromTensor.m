function J = getVectorIndexFromTensor(j,dims)

% dimension of total system
d = prod(dims) ;
n = length(dims) ;

% total index is
J = 1 ;
for k = 1:n
   J = J + (j(k)-1)*prod(dims((k+1):n)) ;
end


end