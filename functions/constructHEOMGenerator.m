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
        n_debye = numel(heom_bath_info.lambda_Ds) ;
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
        cs_array_debye = [] ;
        nus_array_debye = [] ;
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

    if (numel(heom_bath_info.lambda_Ds_pade)>0)
        n_debye_pade_baths = numel(heom_bath_info.lambda_Ds_pade) ;
        mode_info.N_pade = [] ;
        nus_pade = {} ;
        cs_pade = {} ;
        cbars_pade = {} ;
        for j = 1:n_debye_pade_baths 
            [nus_array,cs_array,cbars_array,Delta] = constructPadeDecomp(heom_bath_info.omega_Ds_pade(j),...
                heom_bath_info.lambda_Ds_pade(j),heom_bath_info.beta,heom_bath_info.N_pade(j),heom_bath_info.pade_approximants(j)) ;
            nus_pade{j} = nus_array ; cs_pade{j} = cs_array ; cbars_pade{j} = cs_array ;
            nus = [nus,nus_array] ;
            cs = [cs,cs_array] ;
            cbars = [cbars,cs_array] ;
            Delta_pade(j) = Delta ;
            lambdas = [lambdas, heom_bath_info.lambda_Ds_pade(j)] ;
            mode_info.N_pade = heom_bath_info.N_pade(j) ; 
        end
        n_debye_pade = sum((mode_info.N_pade+1)) ;
        
       
        mode_info.n_debye_pade = n_debye_pade ;
        

    else
        mode_info.n_debye_pade = 0 ;
        mode_info.N_pade = [] ;
        cs_array_debye_pade = [] ;
        nus_array_debye_pade = [] ;
    end

    if (heom_truncation_info.truncation_method == "frequency cut-off")
        % construct the hierarchy structure with frequency truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes,ado_indices_term,modes_term,term_indices] = generateHierarchyFreqTrunc(heom_truncation_info.Gamma_cut,nus);
    elseif (heom_truncation_info.truncation_method == "depth cut-off")
        % construct the hierarchy structure with depth (L) based truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes,ado_indices_term,modes_term,term_indices] = generateHierarchyDepthTrunc(L_max,nus) ;
    elseif (heom_truncation_info.truncation_method == "coupling weighted cut-off")
        % construct the hierarchy structure with coupling weighted depth
        % truncation
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes,ado_indices_term,modes_term,term_indices] = generateHierarchyCouplingWeightedCutoffTrunc(heom_truncation_info.L_cut,heom_truncation_info.p,nus,cs) ;
    elseif (heom_truncation_info.truncation_method == "lambda weighted cut-off")
        % uses the lambda weighted cut-off
        [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes,ado_indices_term,modes_term,term_indices] = generateHierarchyLambdaWeightedCutoffTrunc(heom_truncation_info.L_cut,heom_truncation_info.p,nus,cs,sqrt(lambdas/heom_bath_info.beta));
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
    %     heom_structure.cs_array_debye = cs_array_debye ;
    %     heom_structure.nus_array_debye = nus_array_debye ;
    heom_structure.mode_info = mode_info ;
    heom_structure.ado_indices_term = ado_indices_term ;
    heom_structure.modes_term = modes_term ;
    heom_structure.term_indices = term_indices;
    heom_structure.cs_array_debye = cs_array_debye ;
    heom_structure.nus_array_debye = nus_array_debye ;
    heom_structure.mode_info = mode_info ;
    heom_structure.ado_indices_term = ado_indices_term ;
    heom_structure.modes_term = modes_term ;
    heom_structure.term_indices = term_indices;
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
    ado_indices_term = heom_structure.ado_indices_term ;
    modes_term = heom_structure.modes_term ;
    term_indices = heom_structure.term_indices;
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
n_debye_pade_baths = numel(heom_bath_info.lambda_Ds_pade) ;
n_couplings = size(lower_indices,1) ;

fprintf('N_ado = %d, M = %d\n',[n_ados,M]) ;

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
Xi = sparse([],[],[],d_heom,d_heom) ;
if (heom_truncation_info.heom_termination == "markovian"  || heom_truncation_info.heom_termination == "low temp correction")
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
    Xi = kron(id_ados,Xi) ;
elseif (heom_truncation_info.heom_termination == "NZ2" || heom_truncation_info.heom_termination == "low temp correction NZ2"...
        || heom_truncation_info.heom_termination == "RF2" || heom_truncation_info.heom_termination == "low temp correction RF2" ...
        || heom_truncation_info.heom_termination == "partial resummed")
    k_max = heom_truncation_info.termination_k_max ;
    % get the eigenvalues and eigenvectors of the system liouvillian
    [Pi_sys,lambda_sys] = eig(full(L_sys),'vector') ;
    %     lambda_sys = 0*lambda_sys ;
    Pi_sys_inv = inv(Pi_sys) ;
    cs_term = {} ;
    cbars_term = {} ;
    nus_term = {} ;
    V_comm_Pi_sys = {} ;
    Pi_sys_inv_V_R = {} ;
    Pi_sys_inv_V_L = {} ;
    Pi_sys_inv_V_R_Pi_sys = {} ;
    Pi_sys_inv_V_L_Pi_sys = {} ;
    n_modes_term = [] ;
    if (numel(heom_bath_info.lambda_Ds)>0)
        [nus_array_debye_term,cs_array_debye_term] = generateNusAndCsDebye(heom_bath_info.omega_Ds,...
            heom_bath_info.lambda_Ds,heom_bath_info.beta,k_max) ;
        n_debye_term = (k_max-M)*n_debye_baths ;
        mode_info.n_debye_term = n_debye_term ;
        for j = 1:n_debye_baths
            j_bath = j;
            nus_term{j_bath} = reshape(transpose(nus_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
            cs_term{j_bath} = reshape(transpose(cs_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
            cbars_term{j_bath} = reshape(transpose(cs_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
        end
    else
        mode_info.n_debye_term = 0 ;
    end
    if (numel(heom_bath_info.lambda_OBOs)>0)
        [nus_array_OBO_term,cs_array_OBO_term] = generateNusAndCsOBO(heom_bath_info.gamma_OBOs,...
            heom_bath_info.lambda_OBOs,heom_bath_info.Omega_OBOs,heom_bath_info.beta,k_max) ;
        n_OBO_term = (k_max-M)*n_OBO_baths ;
        mode_info.n_obo_term = n_OBO_term ;
        for j = 1:n_OBO_baths
            j_bath = j+n_debye_baths ;
            nus_term{j_bath} = reshape(transpose(nus_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
            cs_term{j_bath} = reshape(transpose(cs_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
            cbars_term{j_bath} = reshape(transpose(cs_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
        end


    else
        mode_info.n_obo_term = 0 ;
    end
    if (numel(heom_bath_info.lambda_UBOs)>0)
        [nus_array_UBO_term,cs_array_UBO_term,cbars_array_UBO_term] = generateNusAndCsUBO(heom_bath_info.gamma_UBOs,...
            heom_bath_info.lambda_UBOs,heom_bath_info.Omega_UBOs,heom_bath_info.beta,k_max) ;
        n_UBO_term = (k_max-M)*n_UBO_baths ;
        mode_info.n_ubo_term = n_UBO_term ;
        
        for j = 1:n_UBO_baths
            j_bath = j+n_debye_baths+n_OBO_baths ;
            nus_term{j_bath} = reshape(transpose(nus_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
            cs_term{j_bath} = reshape(transpose(cs_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
            cbars_term{j_bath} = reshape(transpose(cbars_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
            
        end
        
    else
        mode_info.n_ubo_term = 0 ;
    end
    if (numel(heom_bath_info.lambda_Ds_pade)>0)
        n_debye_pade_term = 0 ;
        for j = 1:n_debye_pade_baths
            j_bath = j+n_debye_baths+n_OBO_baths+n_UBO_baths ;
            [nus_array_padecorr,cs_array_padecorr] = generateNusAndCsDebye(heom_bath_info.omega_Ds_pade(j),...
                heom_bath_info.lambda_Ds_pade(j),heom_bath_info.beta,k_max) ;
            nus_term{j_bath} = [nus_pade{j}  ,  nus_array_padecorr ] ;
            cs_term{j_bath} = [-cs_pade{j},cs_array_padecorr] ;
            cbars_term{j_bath} = [-cs_pade{j},cs_array_padecorr] ;
            n_debye_pade_term = n_debye_pade_term + numel(nus_term{j_bath}) ;
        end
        mode_info.n_debye_pade_term = n_debye_pade_term ;
    else
        mode_info.n_debye_pade_term = 0 ;
    end


    n_NZ2_term = mode_info.n_debye_term + mode_info.n_obo_term + mode_info.n_ubo_term + mode_info.n_debye_pade_term ;

    for j = 1:n_baths
        V_comm_Pi_sys{j} = V_comm{j} * Pi_sys ;
        Pi_sys_inv_V_R{j} = Pi_sys_inv * V_R{j} ;
        Pi_sys_inv_V_L{j} = Pi_sys_inv * V_L{j} ;
        Pi_sys_inv_V_R_Pi_sys{j} = Pi_sys_inv_V_R{j} * Pi_sys ;
        Pi_sys_inv_V_L_Pi_sys{j} = Pi_sys_inv_V_L{j} * Pi_sys ;
    end
    % create empty Xi
    if (heom_truncation_info.heom_termination == "NZ2" || heom_truncation_info.heom_termination == "low temp correction NZ2"...
            || heom_truncation_info.heom_termination == "partial resummed")
        Xi = sparse([],[],[],d_heom,d_heom,n_ados*d_liou*d_liou) ;
        for J = 1:n_ados
            J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
            for j =1:n_baths
                n_jks_term = zeros(size(cs_term{j})) ;
                Xi(J_block,J_block) = Xi(J_block,J_block) ...
                    - V_comm_Pi_sys{j} * ((sum(((n_jks_term+1).*cs_term{j})./((ado_gammas(J)+nus_term{j})-lambda_sys ),2).* Pi_sys_inv_V_L{j})...
                    - (sum(((n_jks_term+1).*conj(cbars_term{j}))./((ado_gammas(J)+nus_term{j})-lambda_sys ),2).* Pi_sys_inv_V_R{j}));
            end
        end
    elseif (heom_truncation_info.heom_termination == "RF2" || heom_truncation_info.heom_termination == "low temp correction RF2")
        fprintf("adding RF term.\n")
%         Xi_RF = zeros([d_liou,d_liou]) ;
        Xi_RF = sparse([],[],[],d_liou,d_liou) ;
        Delta_lambda_sys = lambda_sys - transpose(lambda_sys) ;
        for j =1:n_baths
            A_j = zeros([d_liou,d_liou]) ;
            Abar_j = zeros([d_liou,d_liou]) ;
            for k = 1:numel(nus_term{j})
                J_jk = 1./(nus_term{j}(k)-  Delta_lambda_sys) ;
                Abar_j = Abar_j + cbars_term{j}(k) * J_jk ;
                A_j = A_j + cs_term{j}(k) * J_jk ;
            end
            Xi_RF = Xi_RF - V_comm_Pi_sys{j} * (...
                (A_j .* Pi_sys_inv_V_L_Pi_sys{j} - Abar_j .* Pi_sys_inv_V_R_Pi_sys{j})) * Pi_sys_inv ;
        end
        Xi = kron(id_ados,Xi_RF) ;
    end
    % use the standard Markovian approximation for all modes with k>k_max ;
%     Xi_markov = zeros([d_liou,d_liou]) ;
    Xi_markov = sparse([],[],[],d_liou,d_liou) ;
    for j = 1:n_debye_baths
        R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) ...
            - sum([cs_array_debye(j,2:end),cs_term{j}]./[nus_array_debye(j,2:end),nus_term{j}]) ;
        Xi_markov = Xi_markov - R_j * V_comm{j}*V_comm{j} ;
    end
    for j = 1:n_OBO_baths
        j_OBO = j + n_debye_baths ;
        ks = (k_max+1):1:max([20*M,100]) ;
        R_j = sum(calculateCkBOs(heom_bath_info.gamma_OBOs(j),heom_bath_info.Omega_OBOs(j),beta,heom_bath_info.lambda_OBOs(j),ks)) ;
        Xi_markov = Xi_markov - R_j * V_comm{j_OBO}*V_comm{j_OBO} ;
    end
    for j = 1:n_UBO_baths
        j_UBO = j + n_debye_baths + n_OBO_baths ;
        ks = (k_max+1):1:max([20*M,100]) ;
        R_j = sum(calculateCkBOs(heom_bath_info.gamma_UBOs(j),heom_bath_info.Omega_UBOs(j),beta,heom_bath_info.lambda_UBOs(j),ks)) ;
        Xi_markov = Xi_markov - R_j * V_comm{j_UBO}*V_comm{j_UBO} ;
    end
    Xi = Xi + kron(id_ados,Xi_markov) ;

end
% add the Pade white-noise term
% use the standard Markovian approximation for all modes with k>k_max ;
if ((numel(heom_bath_info.lambda_Ds_pade)>0)&&(~(heom_truncation_info.heom_termination == "RF2") ...
        && ~(heom_truncation_info.heom_termination == "NZ2") && ~(heom_truncation_info.heom_termination == "low temp correction RF2")...
        && ~(heom_truncation_info.heom_termination == "low temp correction NZ2") ))
    fprintf('Adding white noise Pade term.\n')
%     Xi_pade = zeros([d_liou,d_liou]) ;
    Xi_pade = sparse([],[],[],d_liou,d_liou) ;
    for j = 1:n_debye_pade_baths      
        j_debye_pade = j + n_debye_baths + n_OBO_baths + n_UBO_baths ;
        R_j = Delta_pade(j) ;                  
        Xi_pade = Xi_pade - R_j * V_comm{j_debye_pade}*V_comm{j_debye_pade} ;
    end
    Xi = Xi + kron(id_ados,Xi_pade) ;

end

% add the free system evolution term to the HEOM generator
L_heom = kron(id_ados,L_sys)+Xi ;

% add the decay terms -sum_jk n_jk nu_jk
L_heom = L_heom - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou) ;

heom_structure.Xi = Xi ;
heom_structure.L_sys = L_sys ;
heom_structure.V = V ;
heom_structure.V_comm = V_comm ;
heom_structure.V_L = V_L ;
heom_structure.V_R = V_R ;
heom_structure.d_heom = d_heom ;
heom_structure.n_ados = n_ados ;
heom_structure.d_liou = d_liou ;

% construct the HEOM generator
for r = 1:n_couplings
    J = lower_indices(r) ;
    K = upper_indices(r) ;
    jk_coup = coupled_mode_indices(r) ;
    % need to fix getting the bath index!!!!
    %     j_coup = ceil(jk_coup/(M+1)) ; % get the bath index that is coupling J & K
    j_coup = coupled_bath_indices(r) ;
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

% ORIGINAL INCORRECT TERMINATOR CODE - diagonal only
diag_only_term = false ;
if (isfield(heom_truncation_info,'diagonal_only_term'))
    if (heom_truncation_info.diagonal_only_term)
        diag_only_term = true ;
    end
end

if (heom_truncation_info.heom_termination == "markovian" && diag_only_term)
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
        L_heom(J_block,J_block) = L_heom(J_block,J_block) ...
            -(1.0i/(ado_gammas(J)+nus(jk_term))) * V_comm{j_term}*...
            ((-1.0i*c_jk_term*(n_jk_term+1)) * V_L{j_term} ...
            + (1.0i*conj(cbar_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
        %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/sum(nus.*n_jks_term)) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
        %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;


        %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -...
        %         (1.0i) * V_comm{j_term}*((-L_sys-Xi+(ado_gammas(J)+nus(jk_term))*id_liou)\((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
        %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} )) ;
    end

elseif (heom_truncation_info.heom_termination == "NZ2" && diag_only_term)
    % add a perturbative correction for the ados at which the hierarchy is terminated
    fprintf('Adding NZ2 terminators.\n')
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
        nu_jk_term = nus(jk_term) ;
        J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
        n_jks_term = n_jks ;
        n_jks_term(jk_term) = n_jks_term(jk_term) + 1 ;
        L_heom(J_block,J_block) = L_heom(J_block,J_block) ...
            - V_comm_Pi_sys{j_term} * (((((n_jk_term+1).*c_jk_term)./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_L{j_term})...
            - ((((n_jk_term+1).*conj(cbar_jk_term))./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_R{j_term}));
    end
elseif (heom_truncation_info.heom_termination == "RF2" && diag_only_term)
    % add a perturbative correction for the ados at which the hierarchy is terminated
    fprintf('Adding RF2 terminators.\n')
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
        nu_jk_term = nus(jk_term) ;
        J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;


        J_jk = (n_jk_term+1)./(nu_jk_term-  Delta_lambda_sys) ;
        Abar_j = cbar_jk_term * J_jk ;
        A_j =  c_jk_term * J_jk ;
        L_heom(J_block,J_block) = L_heom(J_block,J_block) ...
             - V_comm_Pi_sys{j_term} * (...
                (A_j .* Pi_sys_inv_V_L_Pi_sys{j_term} - Abar_j .* Pi_sys_inv_V_R_Pi_sys{j_term})) * Pi_sys_inv ;
    end
elseif (heom_truncation_info.heom_termination == "partial resummed")
    % add a perturbative correction for the ados at which the hierarchy is terminated
    [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
    terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
    n_term = size(terminator_ado_indices,1) ;
    n_max = heom_truncation_info.n_max_resum ;
    for r = 1:n_term
        jk_term = terminator_modes(r) ;
        j_term = terminator_bath_indices(r) ;
        J = terminator_ado_indices(r) ;
        n_jks = ado_indices(J,:) ;
        n_jk_term = n_jks(jk_term) ;
        c_jk_term = cs(jk_term) ;
        cbar_jk_term = cbars(jk_term) ;
        nu_jk_term = nus(jk_term) ;
        J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
        T = computeResummedTerminator(V_L{j_term},V_R{j_term},c_jk_term,cbar_jk_term,nu_jk_term,n_jk_term,ado_gammas(J),n_max,Pi_sys,Pi_sys_inv,lambda_sys,L_sys) ;
        L_heom(J_block,J_block) = L_heom(J_block,J_block) + T ;
%         T_NZ = - V_comm_Pi_sys{j_term} * (((((n_jk_term+1).*c_jk_term)./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_L{j_term})...
%             - ((((n_jk_term+1).*conj(cbar_jk_term))./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_R{j_term}))
        
    end
end

if ((heom_truncation_info.heom_termination == "markovian" || heom_truncation_info.heom_termination == "NZ2")&& ~diag_only_term)
    % add a perturbative correction for the ados at which the hierarchy is terminated
    %     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
    %     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
    %     n_term = size(terminator_ado_indices,1) ;
    % get the number of terminating ADOs
    n_term_ados = size(ado_indices_term,1) ;
    for r = 1:n_term_ados
        % get the number of connections up into the explicit hierarchy
        n_term_indices = numel(term_indices{r}) ;
        for ind_J = 1:n_term_indices
            ind_K_range = 1:n_term_indices ;
            %             if (isfield(heom_truncation_info,'diagonal_only_term'))
            %                 if (heom_truncation_info.diagonal_only_term)
            %                     ind_K_range = ind_J ;
            %                 end
            %             end

            for ind_K = ind_K_range
                % get the hierarchy indices of the coupled terms in the
                % truncated hierarchy
                J = term_indices{r}(ind_J) ;
                K = term_indices{r}(ind_K) ;
                % get the block indices
                J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
                K_block = ((K-1)*(d_liou)+1):(K*d_liou) ;
                % get the terminating modes
                term_mode_J = modes_term{r}(ind_J) ;
                term_mode_K = modes_term{r}(ind_K) ;
                % get terminating bath indices
                j_J = getBathIndexFromModeIndex(term_mode_J,mode_info) ;
                j_K = getBathIndexFromModeIndex(term_mode_K,mode_info) ;
                % the the coupling coefficients
                c_J = cs(term_mode_J) ;
                c_K = cs(term_mode_K) ;
                cbar_K = cbars(term_mode_K) ;
                % get the mode excitation levels
                n_J = ado_indices(J,term_mode_J) ;
                n_K = ado_indices(K,term_mode_K) ;
                gamma_ado_term = sum(nus.*ado_indices_term(r,:)) ;
                if (heom_truncation_info.heom_termination == "markovian")
                    L_heom(J_block,K_block) = L_heom(J_block,K_block) ...
                        -(sqrt((n_J+1)*abs(c_J))/gamma_ado_term) * V_comm{j_J}*...
                        (c_K*sqrt((n_K+1)/abs(c_K)) * V_L{j_K} ...
                        - conj(cbar_K)*sqrt((n_K+1)/abs(c_K)) * V_R{j_K} ) ;
                elseif (heom_truncation_info.heom_termination == "NZ2")
                    L_heom(J_block,K_block) = L_heom(J_block,K_block) ...
                        -(sqrt((n_J+1)*abs(c_J))) * V_comm_Pi_sys{j_J}*...
                        ( (c_K*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)).* Pi_sys_inv_V_L{j_K} ...
                        - (conj(cbar_K)*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)) .* Pi_sys_inv_V_R{j_K} ) ;
                end
            end
        end
    end

end



end