function [L_heom_AB] = constructHEOMABGenerator(H_sys_A,H_sys_B,V_As,V_Bs,heom_bath_info,heom_truncation_info,heom_structure_in)


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

% get some dimensions of things for setting up the HEOM generator
d_sys_A = size(H_sys_A,1) ;
d_sys_B = size(H_sys_B,1) ;
d_liou = d_sys_A*d_sys_B ;
n_ados = size(ado_indices,1) ;
d_heom = d_liou*n_ados ; % total dimension of the HEOM system of ADOs
beta = heom_bath_info.beta ;
n_baths = length(V_As) ;
n_debye_baths = numel(heom_bath_info.lambda_Ds) ;
n_OBO_baths = numel(heom_bath_info.lambda_OBOs) ;
n_UBO_baths = numel(heom_bath_info.lambda_UBOs) ;
n_couplings = size(lower_indices,1) ;

% first make some useful system superoperators
id_sys_A = speye(d_sys_A) ;
id_sys_B = speye(d_sys_B) ;
id_liou = speye(d_liou) ;
id_ados = speye(n_ados) ;

% free system liouvillian
L_sys_AB = -1.0i*(kron(H_sys_A,id_sys_B) - kron(id_sys_A,transpose(H_sys_B))) ;

% superoperators that correspond to left (L) and right (R) multiplication
% by the bath operators and the renormalisation operator
V_L = {} ; V_R = {} ; V_comm = {} ;
for j = 1:n_baths
    V_L{j} = kron(V_As{j},id_sys_B) ;
    V_R{j} = kron(id_sys_A,transpose(V_Bs{j})) ;
    V_comm{j} = V_L{j}-V_R{j} ;
end

% % add matsurbara truncation correction
% Xi = sparse([],[],[],d_liou,d_liou) ;
% if (heom_truncation_info.heom_termination == "markovian" )
%     for j = 1:n_debye_baths
%         R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) - sum(cs_array_debye(j,2:end)./nus_array_debye(j,2:end)) ;
%         Xi = Xi - R_j * V_comm{j}*V_comm{j} ;
%     end
%     for j = 1:n_OBO_baths
%         j_OBO = j + n_debye_baths ;
%         ks = (M+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_OBOs(j),heom_bath_info.Omega_OBOs(j),beta,heom_bath_info.lambda_OBOs(j),ks)) ;
%         Xi = Xi - R_j * V_comm{j_OBO}*V_comm{j_OBO} ;
%     end
%     for j = 1:n_UBO_baths
%         j_UBO = j + n_debye_baths + n_OBO_baths ;
%         ks = (M+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_UBOs(j),heom_bath_info.Omega_UBOs(j),beta,heom_bath_info.lambda_UBOs(j),ks)) ;
%         Xi = Xi - R_j * V_comm{j_UBO}*V_comm{j_UBO} ;
%     end
% end
% % add the free system evolution term to the HEOM generator
% L_heom_AB = kron(id_ados,L_sys_AB+Xi) ;

% add matsurbara truncation correction
% Xi = sparse([],[],[],d_liou,d_liou) ;
% if (heom_truncation_info.heom_termination == "markovian" )
%     for j = 1:n_debye_baths
%         R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) - sum(cs_array_debye(j,2:end)./nus_array_debye(j,2:end)) ;
%         Xi = Xi - R_j * V_comm{j}*V_comm{j} ;
%     end
%     for j = 1:n_OBO_baths
%         j_OBO = j + n_debye_baths ;
%         ks = (M+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_OBOs(j),heom_bath_info.Omega_OBOs(j),beta,heom_bath_info.lambda_OBOs(j),ks)) ;
%         Xi = Xi - R_j * V_comm{j_OBO}*V_comm{j_OBO} ;
%     end
%     for j = 1:n_UBO_baths
%         j_UBO = j + n_debye_baths + n_OBO_baths ;
%         ks = (M+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_UBOs(j),heom_bath_info.Omega_UBOs(j),beta,heom_bath_info.lambda_UBOs(j),ks)) ;
%         Xi = Xi - R_j * V_comm{j_UBO}*V_comm{j_UBO} ;
%     end
%     Xi = kron(id_ados,Xi) ;
% elseif (heom_truncation_info.heom_termination == "NZ2")
%     k_max = heom_truncation_info.termination_k_max ;
%     % get the eigenvalues and eigenvectors of the system liouvillian
%     [Pi_sys,lambda_sys] = eig(full(L_sys_AB),'vector') ;
% %     lambda_sys = 0*lambda_sys ;
%     Pi_sys_inv = inv(Pi_sys) ;
%     cs_term = {} ;
%     cbars_term = {} ;
%     nus_term = {} ;
%     V_comm_Pi_sys = {} ;
%     Pi_sys_inv_V_R = {} ;
%     Pi_sys_inv_V_L = {} ;
%     if (numel(heom_bath_info.lambda_Ds)>0)
%         [nus_array_debye_term,cs_array_debye_term] = generateNusAndCsDebye(heom_bath_info.omega_Ds,...
%             heom_bath_info.lambda_Ds,heom_bath_info.beta,k_max) ;
%         n_debye_term = (k_max-M)*n_debye_baths ;
%         mode_info.n_debye_term = n_debye_term ;
%         for j = 1:n_debye_baths
%             j_bath = j;
%             nus_term{j_bath} = reshape(transpose(nus_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
%             cs_term{j_bath} = reshape(transpose(cs_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
%             cbars_term{j_bath} = reshape(transpose(cs_array_debye_term(j,(M+2):end)),[1,(k_max-M)]) ;
%         end
%     else
%         mode_info.n_debye_term = 0 ;
%     end
%     if (numel(heom_bath_info.lambda_OBOs)>0)
%         [nus_array_OBO_term,cs_array_OBO_term] = generateNusAndCsOBO(heom_bath_info.gamma_OBOs,...
%             heom_bath_info.lambda_OBOs,heom_bath_info.Omega_OBOs,heom_bath_info.beta,k_max) ;
%         n_OBO_term = (k_max-M)*n_OBO_baths ;
%         mode_info.n_obo_term = n_OBO_term ;
%         for j = 1:n_OBO_baths
%             j_bath = j+n_debye_baths ;
%             nus_term{j_bath} = reshape(transpose(nus_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%             cs_term{j_bath} = reshape(transpose(cs_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%             cbars_term{j_bath} = reshape(transpose(cs_array_OBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%         end
% 
% 
%     else
%         mode_info.n_obo_term = 0 ;
%     end
%     if (numel(heom_bath_info.lambda_UBOs)>0)
%         [nus_array_UBO_term,cs_array_UBO_term,cbars_array_UBO_term] = generateNusAndCsUBO(heom_bath_info.gamma_UBOs,...
%             heom_bath_info.lambda_UBOs,heom_bath_info.Omega_UBOs,heom_bath_info.beta,k_max) ;
%         n_UBO_term = (k_max-M)*n_UBO_baths ;
%         mode_info.n_ubo_term = n_UBO_term ;
%         for j = 1:n_UBO_baths
%             j_bath = j+n_debye_baths+n_OBO_baths ;
%             nus_term{j_bath} = reshape(transpose(nus_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%             cs_term{j_bath} = reshape(transpose(cs_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%             cbars_term{j_bath} = reshape(transpose(cbars_array_UBO_term(j,(M+3):end)),[1,(k_max-M)]) ;
%         end
%     else
%         mode_info.n_ubo_term = 0 ;
%     end
%     n_NZ2_term = mode_info.n_debye_term + mode_info.n_obo_term + mode_info.n_ubo_term ;
% 
%     for j = 1:n_baths
%         V_comm_Pi_sys{j} = V_comm{j} * Pi_sys ;
%         Pi_sys_inv_V_R{j} = Pi_sys_inv * V_R{j} ;
%         Pi_sys_inv_V_L{j} = Pi_sys_inv * V_L{j} ;
%     end
%     % create empty Xi
%     Xi = sparse([],[],[],d_heom,d_heom,n_ados*d_liou*d_liou) ;
%     for J = 1:n_ados
%         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%         for j =1:n_baths
%             n_jks_term = zeros([1,k_max-M]) ;
%             Xi(J_block,J_block) = Xi(J_block,J_block) ...
%                 - V_comm_Pi_sys{j} * ((sum(((n_jks_term+1).*cs_term{j})./((ado_gammas(J)+nus_term{j})-lambda_sys ),2).* Pi_sys_inv_V_L{j})...
%                 - (sum(((n_jks_term+1).*conj(cbars_term{j}))./((ado_gammas(J)+nus_term{j})-lambda_sys ),2).* Pi_sys_inv_V_R{j}));
%         end
%     end
%     % use the standard Markovian approximation for all modes with k>k_max ;
%     Xi_markov = zeros([d_liou,d_liou]) ;
%     for j = 1:n_debye_baths
%         R_j = 2.0*heom_bath_info.lambda_Ds(j)/(beta*heom_bath_info.omega_Ds(j)) - heom_bath_info.lambda_Ds(j)*cot(beta*heom_bath_info.omega_Ds(j)/2) ...
%             - sum([cs_array_debye(j,2:end),cs_term{j}]./[nus_array_debye(j,2:end),nus_term{j}]) ;
%         Xi_markov = Xi_markov - R_j * V_comm{j}*V_comm{j} ;
%     end
%     for j = 1:n_OBO_baths
%         j_OBO = j + n_debye_baths ;
%         ks = (k_max+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_OBOs(j),heom_bath_info.Omega_OBOs(j),beta,heom_bath_info.lambda_OBOs(j),ks)) ;
%         Xi_markov = Xi_markov - R_j * V_comm{j_OBO}*V_comm{j_OBO} ;
%     end
%     for j = 1:n_UBO_baths
%         j_UBO = j + n_debye_baths + n_OBO_baths ;
%         ks = (k_max+1):1:max([20*M,100]) ;
%         R_j = sum(calculateCkBOs(heom_bath_info.gamma_UBOs(j),heom_bath_info.Omega_UBOs(j),beta,heom_bath_info.lambda_UBOs(j),ks)) ;
%         Xi_markov = Xi_markov - R_j * V_comm{j_UBO}*V_comm{j_UBO} ;
%     end
%     Xi = Xi + kron(id_ados,Xi_markov) ;
%     
% end
% % add the free system evolution term to the HEOM generator
% L_heom_AB = kron(id_ados,L_sys_AB)+Xi ;
L_heom_AB = kron(id_ados,L_sys_AB);
% L_heom_AB = L_heom_AB+Xi ;

% add the decay terms -sum_jk n_jk nu_jk
L_heom_AB = L_heom_AB - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou) ;

% heom_structure.Xi = Xi ;
% heom_structure.L_sys = L_sys_AB ;
% heom_structure.V = V ;
% heom_structure.V_comm = V_comm ;
% heom_structure.V_L = V_L ;
% heom_structure.V_R = V_R ;

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
    L_heom_AB(J_block,K_block) = L_heom_AB(J_block,K_block)...
        -1.0i *sqrt(n_jk_coup*abs(c_jk_coup))* V_comm{j_coup} ;
    % add the coupling from the ADO shallower in the heirarchy
    if (abs(c_jk_coup)>0)
        L_heom_AB(K_block,J_block) = L_heom_AB(K_block,J_block) ...
            + (-1.0i*c_jk_coup*sqrt(n_jk_coup/abs(c_jk_coup))) * V_L{j_coup} ...
            + (1.0i*conj(cbar_jk_coup)*sqrt(n_jk_coup/abs(c_jk_coup))) * V_R{j_coup} ;
    end
end

% ORIGINAL INCORRECT TERMINATION
% if (heom_truncation_info.heom_termination == "markovian")
%     % add a perturbative correction for the ados at which the hierarchy is terminated
%     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
%     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
%     n_term = size(terminator_ado_indices,1) ;
%     for r = 1:n_term
%         jk_term = terminator_modes(r) ;
%         j_term = terminator_bath_indices(r) ;
%         J = terminator_ado_indices(r) ;
%         n_jks = ado_indices(J,:) ;
%         n_jk_term = n_jks(jk_term) ;
%         c_jk_term = cs(jk_term) ;
%         cbar_jk_term = cbars(jk_term) ;
%         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%         n_jks_term = n_jks ;
%         n_jks_term(jk_term) = n_jks_term(jk_term) + 1 ;
%         L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) -(1.0i/(ado_gammas(J)+nus(jk_term))) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%             + (1.0i*conj(cbar_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
%         %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/sum(nus.*n_jks_term)) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
% 
% 
%         %     L_heom(J_block,J_block) = L_heom(J_block,J_block) -...
%         %         (1.0i) * V_comm{j_term}*((-L_sys-Xi+(ado_gammas(J)+nus(jk_term))*id_liou)\((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} )) ;
%     end
% end

% if (heom_truncation_info.heom_termination == "markovian" || heom_truncation_info.heom_termination == "NZ2")
%     % add a perturbative correction for the ados at which the hierarchy is terminated
% %     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
% %     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
% %     n_term = size(terminator_ado_indices,1) ;
%     % get the number of terminating ADOs
%     n_term_ados = size(ado_indices_term,1) ;
%     for r = 1:n_term_ados
%         % get the number of connections up into the explicit hierarchy
%         n_term_indices = numel(term_indices{r}) ;
%         for ind_J = 1:n_term_indices
%             for ind_K = 1:n_term_indices
%                 % get the hierarchy indices of the coupled terms in the
%                 % truncated hierarchy
%                 J = term_indices{r}(ind_J) ;
%                 K = term_indices{r}(ind_K) ;
%                 % get the block indices
%                 J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%                 K_block = ((K-1)*(d_liou)+1):(K*d_liou) ;
%                 % get the terminating modes
%                 term_mode_J = modes_term{r}(ind_J) ;
%                 term_mode_K = modes_term{r}(ind_K) ;
%                 % get terminating bath indices
%                 j_J = getBathIndexFromModeIndex(term_mode_J,mode_info) ;
%                 j_K = getBathIndexFromModeIndex(term_mode_K,mode_info) ;
%                 % the the coupling coefficients
%                 c_J = cs(term_mode_J) ;
%                 c_K = cs(term_mode_K) ;
%                 cbar_K = cbars(term_mode_K) ;
%                 % get the mode excitation levels
%                 n_J = ado_indices(J,term_mode_J) ;
%                 n_K = ado_indices(K,term_mode_K) ;
%                 gamma_ado_term = sum(nus.*ado_indices_term(r,:)) ;
%                 if (heom_truncation_info.heom_termination == "markovian")
%                     L_heom_AB(J_block,K_block) = L_heom_AB(J_block,K_block) ...
%                         -(sqrt((n_J+1)*abs(c_J))/gamma_ado_term) * V_comm{j_J}*...
%                         (c_K*sqrt((n_K+1)/abs(c_K)) * V_L{j_K} ...
%                         - conj(cbar_K)*sqrt((n_K+1)/abs(c_K)) * V_R{j_K} ) ;
%                 elseif (heom_truncation_info.heom_termination == "NZ2")
%                     L_heom_AB(J_block,K_block) = L_heom_AB(J_block,K_block) ...
%                     -(sqrt((n_J+1)*abs(c_J))) * V_comm_Pi_sys{j_J}*...
%                     ( (c_K*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)).* Pi_sys_inv_V_L{j_K} ...
%                     - (conj(cbar_K)*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)) .* Pi_sys_inv_V_R{j_K} ) ;
%                 end
%             end
%         end
%     end
%        
% end

% ORIGINAL INCORRECT TERMINATOR CODE - diagonal only
% diag_only_term = false ;
% if (isfield(heom_truncation_info,'diagonal_only_term'))
%     if (heom_truncation_info.diagonal_only_term)
%         diag_only_term = true ;
%     end
% end
% 
% if (heom_truncation_info.heom_termination == "markovian" && diag_only_term)
%     % add a perturbative correction for the ados at which the hierarchy is terminated
%     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
%     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
%     n_term = size(terminator_ado_indices,1) ;
%     for r = 1:n_term
%         jk_term = terminator_modes(r) ;
%         j_term = terminator_bath_indices(r) ;
%         J = terminator_ado_indices(r) ;
%         n_jks = ado_indices(J,:) ;
%         n_jk_term = n_jks(jk_term) ;
%         c_jk_term = cs(jk_term) ;
%         cbar_jk_term = cbars(jk_term) ;
%         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%         n_jks_term = n_jks ;
%         n_jks_term(jk_term) = n_jks_term(jk_term) + 1 ;
%         L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) ...
%             -(1.0i/(ado_gammas(J)+nus(jk_term))) * V_comm{j_term}*...
%             ((-1.0i*c_jk_term*(n_jk_term+1)) * V_L{j_term} ...
%             + (1.0i*conj(cbar_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
%         %     L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) -(1.0i/sum(nus.*n_jks_term)) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
% 
% 
%         %     L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) -...
%         %         (1.0i) * V_comm{j_term}*((-L_sys-Xi+(ado_gammas(J)+nus(jk_term))*id_liou)\((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         %         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} )) ;
%     end
% 
% elseif (heom_truncation_info.heom_termination == "NZ2" && diag_only_term)
%     % add a perturbative correction for the ados at which the hierarchy is terminated
%     fprintf('Adding NZ2 terminators.\n')
%     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
%     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
%     n_term = size(terminator_ado_indices,1) ;
%     for r = 1:n_term
%         jk_term = terminator_modes(r) ;
%         j_term = terminator_bath_indices(r) ;
%         J = terminator_ado_indices(r) ;
%         n_jks = ado_indices(J,:) ;
%         n_jk_term = n_jks(jk_term) ;
%         c_jk_term = cs(jk_term) ;
%         cbar_jk_term = cbars(jk_term) ;
%         nu_jk_term = nus(jk_term) ;
%         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%         n_jks_term = n_jks ;
%         n_jks_term(jk_term) = n_jks_term(jk_term) + 1 ;
%         L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) ...
%             - V_comm_Pi_sys{j_term} * (((((n_jk_term+1).*c_jk_term)./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_L{j_term})...
%             - ((((n_jk_term+1).*conj(cbar_jk_term))./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_R{j_term}));
%     end
% % elseif (heom_truncation_info.heom_termination == "RF2" && diag_only_term)
% %     add a perturbative correction for the ados at which the hierarchy is terminated
% %     fprintf('Adding RF2 terminators.\n')
% %     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
% %     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
% %     n_term = size(terminator_ado_indices,1) ;
% %     for r = 1:n_term
% %         jk_term = terminator_modes(r) ;
% %         j_term = terminator_bath_indices(r) ;
% %         J = terminator_ado_indices(r) ;
% %         n_jks = ado_indices(J,:) ;
% %         n_jk_term = n_jks(jk_term) ;
% %         c_jk_term = cs(jk_term) ;
% %         cbar_jk_term = cbars(jk_term) ;
% %         nu_jk_term = nus(jk_term) ;
% %         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
% % 
% % 
% %         J_jk = (n_jk_term+1)./(nu_jk_term-  Delta_lambda_sys) ;
% %         Abar_j = cbar_jk_term * J_jk ;
% %         A_j =  c_jk_term * J_jk ;
% %         L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) ...
% %              - V_comm_Pi_sys{j_term} * (...
% %                 (A_j .* Pi_sys_inv_V_L_Pi_sys{j_term} - Abar_j .* Pi_sys_inv_V_R_Pi_sys{j_term})) * Pi_sys_inv ;
% %     end
% % elseif (heom_truncation_info.heom_termination == "partial resummed")
% %     % add a perturbative correction for the ados at which the hierarchy is terminated
% %     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
% %     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
% %     n_term = size(terminator_ado_indices,1) ;
% %     n_max = heom_truncation_info.n_max_resum ;
% %     for r = 1:n_term
% %         jk_term = terminator_modes(r) ;
% %         j_term = terminator_bath_indices(r) ;
% %         J = terminator_ado_indices(r) ;
% %         n_jks = ado_indices(J,:) ;
% %         n_jk_term = n_jks(jk_term) ;
% %         c_jk_term = cs(jk_term) ;
% %         cbar_jk_term = cbars(jk_term) ;
% %         nu_jk_term = nus(jk_term) ;
% %         J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
% %         T = computeResummedTerminator(V_L{j_term},V_R{j_term},c_jk_term,cbar_jk_term,nu_jk_term,n_jk_term,ado_gammas(J),n_max,Pi_sys,Pi_sys_inv,lambda_sys,L_sys) ;
% %         L_heom_AB(J_block,J_block) = L_heom_AB(J_block,J_block) + T ;
% % %         T_NZ = - V_comm_Pi_sys{j_term} * (((((n_jk_term+1).*c_jk_term)./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_L{j_term})...
% % %             - ((((n_jk_term+1).*conj(cbar_jk_term))./((ado_gammas(J)+nu_jk_term)-lambda_sys )).* Pi_sys_inv_V_R{j_term}))
% %         
% %     end
% end
% 
% if ((heom_truncation_info.heom_termination == "markovian" || heom_truncation_info.heom_termination == "NZ2")&& ~diag_only_term)
%     % add a perturbative correction for the ados at which the hierarchy is terminated
%     %     [terminator_ado_indices,terminator_modes] = find(truncated_coupled_modes) ;
%     %     terminator_bath_indices = getCoupledBathIndices(terminator_modes,mode_info) ;
%     %     n_term = size(terminator_ado_indices,1) ;
%     % get the number of terminating ADOs
%     n_term_ados = size(ado_indices_term,1) ;
%     for r = 1:n_term_ados
%         % get the number of connections up into the explicit hierarchy
%         n_term_indices = numel(term_indices{r}) ;
%         for ind_J = 1:n_term_indices
%             ind_K_range = 1:n_term_indices ;
%             %             if (isfield(heom_truncation_info,'diagonal_only_term'))
%             %                 if (heom_truncation_info.diagonal_only_term)
%             %                     ind_K_range = ind_J ;
%             %                 end
%             %             end
% 
%             for ind_K = ind_K_range
%                 % get the hierarchy indices of the coupled terms in the
%                 % truncated hierarchy
%                 J = term_indices{r}(ind_J) ;
%                 K = term_indices{r}(ind_K) ;
%                 % get the block indices
%                 J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%                 K_block = ((K-1)*(d_liou)+1):(K*d_liou) ;
%                 % get the terminating modes
%                 term_mode_J = modes_term{r}(ind_J) ;
%                 term_mode_K = modes_term{r}(ind_K) ;
%                 % get terminating bath indices
%                 j_J = getBathIndexFromModeIndex(term_mode_J,mode_info) ;
%                 j_K = getBathIndexFromModeIndex(term_mode_K,mode_info) ;
%                 % the the coupling coefficients
%                 c_J = cs(term_mode_J) ;
%                 c_K = cs(term_mode_K) ;
%                 cbar_K = cbars(term_mode_K) ;
%                 % get the mode excitation levels
%                 n_J = ado_indices(J,term_mode_J) ;
%                 n_K = ado_indices(K,term_mode_K) ;
%                 gamma_ado_term = sum(nus.*ado_indices_term(r,:)) ;
%                 if (heom_truncation_info.heom_termination == "markovian")
%                     L_heom_AB(J_block,K_block) = L_heom_AB(J_block,K_block) ...
%                         -(sqrt((n_J+1)*abs(c_J))/gamma_ado_term) * V_comm{j_J}*...
%                         (c_K*sqrt((n_K+1)/abs(c_K)) * V_L{j_K} ...
%                         - conj(cbar_K)*sqrt((n_K+1)/abs(c_K)) * V_R{j_K} ) ;
%                 elseif (heom_truncation_info.heom_termination == "NZ2")
%                     L_heom_AB(J_block,K_block) = L_heom_AB(J_block,K_block) ...
%                         -(sqrt((n_J+1)*abs(c_J))) * V_comm_Pi_sys{j_J}*...
%                         ( (c_K*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)).* Pi_sys_inv_V_L{j_K} ...
%                         - (conj(cbar_K)*sqrt((n_K+1)/abs(c_K))./ (gamma_ado_term-lambda_sys)) .* Pi_sys_inv_V_R{j_K} ) ;
%                 end
%             end
%         end
%     end
% 
% end

end