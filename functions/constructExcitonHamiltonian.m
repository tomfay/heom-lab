function H = constructExcitonHamiltonian(E_sites, V_sites, coupled_indices)
    % constructs an exciton Hamiltonian
    n_sites = length(E_sites) ;
    
    % site energy part
    H = diag(E_sites) ;
    
    % site couplings
    n_coup = length(V_sites) ;
    for i = 1:n_coup
        n = coupled_indices(i,1) ;
        m = coupled_indices(i,2) ;
        H(n,m) = V_sites(i) ;
        H(m,n) = V_sites(i) ;
    end

end