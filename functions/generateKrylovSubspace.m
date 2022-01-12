function [A_krylov, krylov_basis] = generateKrylovSubspace(A, v, dim_krylov)
    % A_krylov is the krylov representation of A 
    % krylov_basis is a (dim_full x dim_krylov) array of the basis vectors
    % of the Krylov subspace. Each vector has a 2-norm equal to 1.

    % obtain the full space dimension and set up the matrix of basis
    % vectors and the matrix representation of A in this space (upper
    % hessenberg form)
    % This is done using the Arnoldi algorithm
    dim_full = size(v,1) ;
    A_krylov = (1.0+0.0i)*zeros(dim_krylov, dim_krylov) ;
    krylov_basis = (1.0+0.0i)*zeros(dim_full,dim_krylov) ;
    
    % first vector is normalised starting vector
    krylov_basis(:,1) = (1.0/norm(v)) *  v ;
    
    tol = 1e-14 ;
    % build the krylov subspace
    for k = 1:(dim_krylov-1)
        % generate new vector
        u = A * krylov_basis(:,k) ;
        % calculate A matrix elements
        A_krylov(1:k,k) = krylov_basis(:,1:k)' * u ;
        % gram-schmidt orthogonalisation
        u = u - krylov_basis(:,1:k) * A_krylov(1:k,k) ;
        % add final element
        A_krylov(k+1,k) = norm(u) ;
        % add the next vector to the basis
        if abs(A_krylov(k+1,k)) > tol
           krylov_basis(:,k+1) = (1.0 / A_krylov(k+1,k)) * u ;
        else
           krylov_basis(:,k+1) = zeros(dim_full,1) ;
        end
        
    end
    % repeat for final vector
    k = dim_krylov ;
    u = A * krylov_basis(:,k) ;
    A_krylov(1:k,k) = krylov_basis(:,1:k)' * u ;
%     u = u - krylov_basis(:,1:k) * A_krylov(1:k,k) ;

end