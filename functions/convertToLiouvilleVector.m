function liouville_vector = convertToLiouvilleVector(operator)
    
    % Hilbert space dimension
    n = size(operator,1) ;
    
    % reshapes operator to give lioville vector
    liouville_vector = reshape(transpose(operator), [n^2,1]) ;

end
