function K_incoh = constructIncoherentRateOperator(incoh_processes,n_ados,d_lious) 

d_heom = sum(n_ados.*d_lious) ;
n_blocks = numel(n_ados) ;
K_incoh = sparse([],[],[],d_heom,d_heom) ;
n_incoh_processes = numel(incoh_processes.rates) ;
block_indices = {};
d_heoms = [] ;
d_hilbs = [] ;
for k = 1:n_blocks
   d_heoms(k) = n_ados(k) * d_lious(k) ;
   d_hilbs(k) = sqrt(d_lious(k)) ;
   block_indices{k} = (1+sum(d_heoms(1:end-1))):sum(d_heoms(1:end)) ;
end
for r = 1:n_incoh_processes
    i = incoh_processes.coupled_blocks(r,1) ;
    f = incoh_processes.coupled_blocks(r,2) ;
    Gamma = incoh_processes.Gammas{r} ;
    rate = incoh_processes.rates(r) ;
    K_incoh(block_indices{i},block_indices{i}) = K_incoh(block_indices{i},block_indices{i}) ...
        - (0.5*rate) * kron(speye(n_ados(i)),kron(Gamma*(Gamma'),speye(d_hilbs(i)))+kron(speye(d_hilbs(i)),transpose(Gamma*(Gamma')))) ;
    K_incoh(block_indices{f},block_indices{i}) = K_incoh(block_indices{f},block_indices{i}) ...
        + rate *  kron(spdiags(ones([n_ados(f),1]),[0],n_ados(f),n_ados(i)),kron(Gamma',transpose(Gamma))) ; 
end

end