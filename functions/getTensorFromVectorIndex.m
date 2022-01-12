function j = getTensorFromVectorIndex(J,dims)

% dimension of total system
d = prod(dims) ;
n = length(dims) ;

% total index is
j = zeros([1,n]) ;
J_red = J-1 ;
for k = 1:n
   d = prod(dims((k+1):n)) ;
   j(k) = floor(J_red/d) + 1 ;
   J_red = J_red - d * (j(k)-1) ;
end


end