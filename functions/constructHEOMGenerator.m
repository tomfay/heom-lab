function [L_heom,heom_structure] = constructHEOMGenerator(H_sys,heom_bath_info,heom_truncation_info,heom_structure_in)

if (nargin~=4)

    % set up the hierarchy using the Djkstra frequency cut-off method using
    % matsubara decompositions of the bath correlation functions
    if (heom_truncation_info.truncation_method == "frequency cut-off")
        % get the maximum matsubara mode that contributrs to the hierarchy
        M = findMaxMj(heom_truncation_info.Gamma_cut, heom_bath_info.beta) ;
    elseif (heom_truncation_info.truncation_method == "depth cut-off")
        % get the specified maximum matsurbara mode
        M = heom_truncation_info.M_max ;
        % get the max depth of the hierarchy
        L_max = heom_truncation_info.L_max ;
    elseif (heom_truncation_info.truncation_method == "coupling weighted cut-off")
        % get the maximum value of M for each mode
        M_max_debye = findMaxMjWeightedCutoff(heom_truncation_info.L_cut, heom_bath_info.omega_Ds, heom_bath_info.lambda_Ds, heom_bath_info.beta) ;
        M_max_OBO = findMaxMjWeightedCutoffBO(heom_truncation_info.L_cut, heom_bath_info.Omega_OBOs, heom_bath_info.lambda_OBOs, heom_bath_info.gamma_OBOs, heom_bath_info.beta) ;
        M_max_UBO = findMaxMjWeightedCutoffBO(heom_truncation_info.L_cut, heom_bath_info.Omega_UBOs, heom_bath_info.lambda_UBOs, heom_bath_info.gamma_UBOs, heom_bath_info.beta) ;
        M = max([M_max_debye,M_max_OBO,M_max_UBO]) ;
    elseif (heom_truncation_info.truncation_method == "lambda weighted cut-off")
        lambda_baths = [heom_bath_info.lambda_Ds, heom_bath_info.lambda_OBOs, heom_bath_info.lambda_UBOs] ;
        M = findMaxMj(heom_truncation_info.L_cut*sqrt(lambda_baths./heom_bath_info.beta), heom_bath_info.beta) ;
    end

    % get arrays the ado frequencies (nus) and coupling coefficents (cs)
    % for the debye, OBO and UBO baths
    nus = [] ;
    cs = [] ;
    cbars = [] ; 
    lambdas = [] ;
    mode_info = struct() ;
    mode_info.M = M ;
    
    if (numel(heom_bath_info.lambda_Ds)>0)
        [nus_array_debye,cs_array_debye] = generateNusAndCsDebye(heom_bath_info.omega_Ds,...
            heom_bath_info.lambda_Ds,heom_bath_info.beta,M) ;
        n_debye = numel(nus_array_debye) ;
        nus = [nus,reshape(transpose(nus_array_debye),[1,n_debye])] ;
        cs = [cs,reshape(transpose(cs_array_debye),[1,n_debye])] ;
        cbars = [cbars,reshape(transpose(cs_array_debye),[1,n_debye])] ;
        mode_info.n_debye = n_debye ;
        lambdas = [lambdas, reshape(repmat(heom_bath_info.lambda_Ds,[M+1,1]),[1,n_debye])] ;

    else
        mode_info.n_debye = 0 ;
    end
    if (numel(heom_bath_info.lambda_OBOs)>0)
        [nus_array_OBO,cs_array_OBO] = generateNusAndCsOBO(heom_bath_info.gamma_OBOs,...
            heom_bath_info.lambda_OBOs,heom_bath_info.Omega_OBOs,heom_bath_info.beta,M) ;
        n_OBO = numel(nus_array_OBO) ;
        nus = [nus,reshape(transpose(nus_array_OBO),[1,n_OBO])] ;
        cs = [cs,reshape(transpose(cs_array_OBO),[1,n_OBO])] ;
        cbars = [cbars,reshape(transpose(cs_array_OBO),[1,n_OBO])] ;
        mode_info.n_obo = n_OBO ;
        lambdas = [lambdas, reshape(repmat(heom_bath_info.lambda_OBOs,[M+2,1]),[1,n_OBO])] ;
    else
        mode_info.n_obo = 0 ;
    end
    if (numel(heom_bath_info.lambda_UBOs)>0)
        [nus_array_UBO,cs_array_UBO,cbars_array_UBO] = generateNusAndCsUBO(heom_bath_info.gamma_UBOs,...
            heom_bath_info.lambda_UBOs,heom_bath_info.Omega_UBOs,heom_bath_info.beta,M) ;
        n_UBO = numel(nus_array_UBO) ;
        nus = [nus,reshape(transpose(nus_array_UBO),[1,n_UBO])] ;
        cs = [cs,reshape(transpose(cs_array_UBO),[1,n_UBO])] ;
        cbars = [cbars,reshape(transpose(cbars_array_UBO),[1,n_UBO])] ;
        mode_info.n_ubo = n_UBO ;
        lambdas = [lambdas, reshape(repmat(heom_bath_info.lambda_UBOs,[M+2,1]),[1,n_UBO])] ;
    else
        mode_info.n_ubo = 0 ;
    end
    if (heom_truncation_info.truncation_method == "frequency cut-off")
        % construct the hierarchy structure with frequency truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyFreqTrunc(heom_truncation_info.Gamma_cut,nus);
    elseif (heom_truncation_info.truncation_method == "depth cut-off")
        % construct the hierarchy structure with depth (L) based truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyDepthTrunc(L_max,nus) ;
    elseif (heom_truncation_info.truncation_method == "coupling weighted cut-off")
        % construct the hierarchy structure with coupling weighted depth
        % truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyCouplingWeightedCutoffTrunc(heom_truncation_info.L_cut,heom_truncation_info.p,nus,cs) ;
    elseif (heom_truncation_info.truncation_method == "lambda weighted cut-off")
        % uses the lambda weighted cut-off
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyLambdaWeightedCutoffTrunc(heom_truncation_info.L_cut,heom_truncation_info.p,nus,cs,sqrt(lambdas/heom_bath_info.beta));
    end

    % create an array of the coupled mode indices
    coupled_bath_indices = getCoupledBathIndices(coupled_mode_indices,mode_info) ;

    heom_structure = struct() ;
    heom_structure.ado_gammas = ado_gammas ;
    heom_structure.ado_indices = ado_indices ;
    heom_structure.lower_indices = lower_indices ;
    heom_structure.upper_indices = upper_indices ;
    heom_structure.coupled_mode_indices = coupled_mode_indices ;
    heom_structure.truncated_coupled_modes = truncated_coupled_modes ;
    heom_structure.coupled_bath_indices = coupled_bath_indices ;
    heom_structure.nus = nus ;
    heom_structure.cs = cs ;
    heom_structure.cbars = cbars ;
    heom_structure.M = M ;
    heom_structure.cs_array_debye = cs_array_debye ;
    heom_structure.nus_array_debye = nus_array_debye ;
    heom_structure.mode_info = mode_info ;
elseif (nargin==4)
    heom_structure = heom_structure_in ;
    ado_gammas = heom_structure.ado_gammas  ;
    ado_indices = heom_structure.ado_indices  ;
    lower_indices = heom_structure.lower_indices  ;
    upper_indices = heom_structure.upper_indices  ;
    coupled_mode_indices = heom_structure.coupled_mode_indices   ;
    truncated_coupled_modes = heom_structure.truncated_coupled_modes  ;
    coupled_bath_indices = heom_structure.coupled_bath_indices  ;
    nus = heom_structure.nus  ;
    cs = heom_structure.cs ;
    cbars = heom_structure.cbars  ;
    M = heom_structure.M ;
    cs_array_debye = heom_structure.cs_array_debye ;
    nus_array_debye = heom_structure.nus_array_debye ;
    mode_info = heom_structure.mode_info ;
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
if (heom_truncation_info.heom_termination == "markovian" )
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
end
% add the free system evolution term to the HEOM generator
L_heom = kron(id_ados,L_sys+Xi) ;

% add the decay terms -sum_jk n_jk nu_jk
L_heom = L_heom - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou) ;

heom_structure.Xi = Xi ;
heom_structure.L_sys = L_sys ;
heom_structure.V = V ;
heom_structure.V_comm = V_comm ;
heom_structure.V_L = V_L ;
heom_structure.V_R = V_R ;

% construct the HEOM generator
for r = 1:n_couplings
    J = lower_indices(r) ;
    K = upper_indices(r) ;
    jk_coup = coupled_mode_indices(r) ;
    % need to fix getting the bath index!!!!
    %     j_coup = ceil(jk_coup/(M+1)) ; % get the bath index that is coupling J & K
    j_coup = coupled_bath_indices(r) ;
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
    if (abs(c_jk_coup)>0)
        L_heom(K_block,J_block) = L_heom(K_block,J_block) ...
            + (-1.0i*c_jk_coup*sqrt(n_jk_coup/abs(c_jk_coup))) * V_L{j_coup} ...
            + (1.0i*conj(cbar_jk_coup)*sqrt(n_jk_coup/abs(c_jk_coup))) * V_R{j_coup} ;
    end
end


if (heom_truncation_info.heom_termination == "markovian")
    % add a perturbative correction for the ados at which the hierarchy is terminated
    [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
    terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
    n_term = size(terminator_ado_indices,1) ;
    for r = 1:n_term
        jk_term = terminator_modes(r) ;
        j_term = terminator_bath_indices(r) ;
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