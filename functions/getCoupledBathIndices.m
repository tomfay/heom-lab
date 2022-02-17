function coupled_bath_indices = getCoupledBathIndices(coupled_mode_indices,mode_info)

n_couplings = numel(coupled_mode_indices) ;
coupled_bath_indices = zeros([1,n_couplings]) ;
n_debye = mode_info.n_debye ;
n_ubo = mode_info.n_ubo ;
n_obo = mode_info.n_obo ;
n_bo = n_ubo + n_obo ;
M = mode_info.M ;
debye_index_max = (M+1)*n_debye ;
for n = 1:n_couplings
    jk_index = coupled_mode_indices(n) ;
    if (jk_index <= debye_index_max)
        coupled_bath_indices(n) = ceil(jk_index/(M+1)) ;
    else
        coupled_bath_indices(n) = n_debye + ceil((jk_index - n_debye)/(M+2)) ;
    end
end

end