function [L_heom,ado_indices] = constructHEOMGenerator(H_sys,heom_bath_info,heom_truncation_info)

% set up the hierarchy using the Djkstra frequency cut-off method using
% matsubara decompositions of the bath correlation functions
if (heom_truncation_info.truncation_method == 'frequency cut-off')
    % get the maximum matsubara mode that contributrs to the hierarchy
    M = findMaxMj(heom_truncation_info.Gamma_cut, heom_bath_info.beta) ;
    % get arrays the ado frequencies (nus) and coupling coefficents (cs)
    % for the debye, OBO and UBO baths
    nus = [] ;
    cs = [] ;
    cbars = [] ;
    if (numel(heom_bath_info.lambda_Ds)>0)
        [nus_array_debye,cs_array_debye] = generateNusAndCsDebye(heom_bath_info.omega_Ds,...
            heom_bath_info.lambda_Ds,heom_bath_info.beta,M) ;
        n_debye = numel(nus_array_debye) ;
        nus = [nus,reshape(transpose(nus_array_debye),[1,n_debye])] ;
        cs = [cs,reshape(transpose(cs_array_debye),[1,n_debye])] ;
        cbars = [cbars,reshape(transpose(cs_array_debye),[1,n_debye])] ;
    end
    if (numel(heom_bath_info.lambda_OBOs)>0)
        [nus_array_OBO,cs_array_OBO] = generateNusAndCsOBO(heom_bath_info.gamma_OBOs,...
            heom_bath_info.lambda_OBOs,heom_bath_info.Omega_OBOs,heom_bath_info.beta,M) ;
        n_OBO = numel(nus_array_OBO) ;
        nus = [nus,reshape(transpose(nus_array_OBO),[1,n_OBO])] ;
        cs = [cs,reshape(transpose(cs_array_OBO),[1,n_OBO])] ;
        cbars = [cbars,reshape(transpose(cs_array_OBO),[1,n_OBO])] ;
    end
    if (numel(heom_bath_info.lambda_UBOs)>0)
        [nus_array_UBO,cs_array_UBO,cbars_array_UBO] = generateNusAndCsUBO(heom_bath_info.gamma_UBOs,...
            heom_bath_info.lambda_UBOs,heom_bath_info.Omega_UBOs,heom_bath_info.beta,M) ;
        n_UBO = numel(nus_array_UBO) ;
        nus = [nus,reshape(transpose(nus_array_UBO),[1,n_UBO])] ;
        cs = [cs,reshape(transpose(cs_array_UBO),[1,n_UBO])] ;
        cbars = [cbars,reshape(transpose(cbars_array_UBO),[1,n_UBO])] ;
    end

    % construct the hierarchy structure with frequency truncation
    [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyFreqTrunc(heom_truncation_info.Gamma_cut,nus);

end

% get some dimensions of things for setting up the HEOM generator
d_sys = size(H_sys,1) ;
d_liou = d_sys^2 ;
n_ados = size(ado_indices,1) ;
d_heom = d_liou*n_ados ; % total dimension of the HEOM system of ADOs
V = heom_bath_info.Vs ;
beta = heom_bath_info.beta ;
n_baths = length(V) ;
n_debye_baths = numel(heom_bath_info.lambda_Ds) ;
n_OBO_baths = numel(heom_bath_info.lambda_OBOs) ;
n_UBO_baths = numel(heom_bath_info.lambda_UBOs) ;
n_couplings = size(lower_indices,1) ;

% first make some useful system superoperators
id_sys = speye(d_sys) ;
id_liou = speye(d_liou) ;
id_ados = speye(n_ados) ;

% free system liouvillian
L_sys = -1.0i*(kron(H_sys,id_sys) - kron(id_sys,transpose(H_sys))) ;

% superoperators that correspond to left (L) and right (R) multiplication
% by the bath operators and the renormalisation operator
V_L = {} ; V_R = {} ; V_comm = {} ;
for j = 1:n_baths
    V_L{j} = kron(V{j},id_sys) ;
    V_R{j} = kron(id_sys,transpose(V{j})) ;
    V_comm{j} = V_L{j}-V_R{j} ;
end

% add matsurbara truncation correction
Xi = sparse([],[],[],d_liou,d_liou) ;
for j = 1:n_debye_baths
    R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) - sum(cs_array_debye(j,2:end)./nus_array_debye(j,2:end)) ;
    Xi = Xi - R_j * V_comm{j}*V_comm{j} ;
end
for j = 1:n_OBO_baths
    j_OBO = j + n_debye_baths ;
    ks = (M+1):1:max([20*M,100]) ;
    R_j = sum(calculateCkBOs(heom_bath_info.gamma_OBOs(j),heom_bath_info.Omega_OBOs(j),beta,heom_bath_info.lambda_OBOs(j),ks)) ;
    Xi = Xi - R_j * V_comm{j_OBO}*V_comm{j_OBO} ;
end
for j = 1:n_UBO_baths
    j_UBO = j + n_debye_baths + n_OBO_baths ;
    ks = (M+1):1:max([20*M,100]) ;
    R_j = sum(calculateCkBOs(heom_bath_info.gamma_UBOs(j),heom_bath_info.Omega_UBOs(j),beta,heom_bath_info.lambda_UBOs(j),ks)) ;
    Xi = Xi - R_j * V_comm{j_UBO}*V_comm{j_UBO} ;
end

% add the free system evolution term to the HEOM generator
L_heom = kron(id_ados,L_sys+Xi) ;

% add the decay terms -sum_jk n_jk nu_jk
L_heom = L_heom - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou) ; 

% construct the HEOM generator
for r = 1:n_couplings
    J = lower_indices(r) ;
    K = upper_indices(r) ;
    jk_coup = coupled_mode_indices(r) ;
    j_coup = ceil(jk_coup/(M+2)) ; % get the bath index that is coupling J & K
    k_coup = jk_coup - j_coup*(M+1) ; % Matsubara mode index, k=0 is non-Matsubara term
    J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
    K_block = ((K-1)*(d_liou)+1):(K*d_liou) ;
    n_jk_coup = ado_indices(K,jk_coup) ;
    c_jk_coup = cs(jk_coup) ;
    cbar_jk_coup = cbars(jk_coup) ;
    % add the coupling from the ADO deeper in the heirarchy
    % d/dt rho_J = ... -i [V_j , rho_K]
    L_heom(J_block,K_block) = L_heom(J_block,K_block)...
        -1.0i *sqrt(n_jk_coup*abs(c_jk_coup))* V_comm{j_coup} ;
    % add the coupling from the ADO shallower in the heirarchy
    L_heom(K_block,J_block) = L_heom(K_block,J_block) ...
        + (-1.0i*c_jk_coup*sqrt(n_jk_coup/abs(c_jk_coup))) * V_L{j_coup} ...
        + (1.0i*conj(cbar_jk_coup)*sqrt(n_jk_coup/abs(c_jk_coup))) * V_R{j_coup} ;
end


if (heom_truncation_info.heom_termination == 'markovian')
% add a perturbative correction for the ados at which the hierarchy is terminated 
[terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
n_term = size(terminator_ado_indices,1) ;
for r = 1:n_term
    jk_term = terminator_modes(r) ;
    j_term = ceil(jk_term/(M+2)) ;
    J = terminator_ado_indices(r) ;
    n_jks = ado_indices(J,:) ;
    n_jk_term = n_jks(jk_term) ;
    c_jk_term = cs(jk_term) ;
    cbar_jk_term = cbars(jk_term) ;
    J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
    n_jks_term = n_jks ;
    n_jks_term(jk_term) = n_jks_term(jk_term) + 1 ;
    L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/(ado_gammas(J)+nus(jk_term))) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
        + (1.0i*conj(cbar_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
    %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/sum(nus.*n_jks_term)) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
    %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;


    %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -...
    %         (1.0i) * V_comm{j_term}*((-L_sys-Xi+(ado_gammas(J)+nus(jk_term))*id_liou)\((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
    %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} )) ;
end
end
end