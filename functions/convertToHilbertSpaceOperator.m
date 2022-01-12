function operator = convertToHilbertSpaceOperator(liouville_vector) 

    % dimension of hilbert space
    n = sqrt(length(liouville_vector)) ;
    
    % reshapes to give hilbert space operator
    operator = reshape(liouville_vector, [n,n]) ;

end