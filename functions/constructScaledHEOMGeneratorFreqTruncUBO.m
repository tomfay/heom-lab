function [L_heom,ado_indices] = constructScaledHEOMGeneratorFreqTruncUBO(H_sys,V,gamma_Bs,lambda_Bs,Omega_Bs,beta,Gamma_cut)

% first find the maximum M, highest Matsubara mode to be treated explicitly
% given the value of Gamma_cut & beta
M = findMaxMjOBO(Gamma_cut,beta) ;

% find the set of frequencies and expansion coefficients for each bath
% using the Matsubara expansion up to mode M
[nus_array,cs_array, cbars_array] = generateNusAndCsUBO(gamma_Bs,lambda_Bs,Omega_Bs,beta,M) ;
n_modes = numel(nus_array) ;
nus = reshape(transpose(nus_array),[1,n_modes]) ;
cs = reshape(transpose(cs_array),[1,n_modes]) ;
cbars = reshape(transpose(cbars_array),[1,n_modes]) ;


% construct a list of the ADO excitation numbers given the cut-off
% [ado_indices,ado_gammas] = generateADOIndices(Gamma_cut,nus) ;

% find which ados in the truncated set are coupled to each other, and which
% mode couples each ADO
% [lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = findCoupledADOIndices(ado_indices) ;

% construct the hierarchy with frequency truncation
[ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyFreqTrunc(Gamma_cut,nus);


% get some dimensions of things for setting up the HEOM generator
d_sys = size(H_sys,1) ;
d_liou = d_sys^2 ;
n_ados = size(ado_indices,1) ;
d_heom = d_liou*n_ados ; % total dimension of the HEOM system of ADOs 
n_baths = length(V) ;
n_couplings = size(lower_indices,1) ;

% first make some useful system superoperators
id_sys = speye(d_sys) ;
id_liou = speye(d_liou) ;
id_ados = speye(n_ados) ;
% free system liouvillian
L_sys = -1.0i*(kron(H_sys,id_sys) - kron(id_sys,transpose(H_sys))) ;
% superoperators that correspond to left (L) and right (R) multiplication
% by the bath operators and the renormalisation operator
Xi = sparse([],[],[],d_liou,d_liou) ;
V_L = {} ; V_R = {} ; V_comm = {} ;
for j = 1:n_baths
    V_L{j} = kron(V{j},id_sys) ;
    V_R{j} = kron(id_sys,transpose(V{j})) ;
    V_comm{j} = V_L{j}-V_R{j} ;
%     R_j = 2.0*lambda_Bs(j)/(beta*gamma_Bs(j)) - lambda_Bs(j)*cot(beta*gamma_Bs(j)/2) - sum(cs_array(j,2:end)./nus_array(j,2:end)) ;
    % need to correct the renormalisation for OBO
%     R_j = sparse([],[],[],d_liou,d_liou) ; 
    ks = (M+1):1:max([20*M,100]) ;
    R_j = sum(calculateCkBOs(gamma_Bs(j),Omega_Bs(j),beta,lambda_Bs(j),ks)) ;
    Xi = Xi - R_j * V_comm{j}*V_comm{j} ;
end


% add the free system evolution term to the HEOM generator
L_heom = kron(id_ados,L_sys+Xi) ;

% add the decay terms -sum_jk n_jk nu_jk
L_heom = L_heom - kron(spdiags([ado_gammas],[0],n_ados,n_ados),id_liou) ; 

% add the coupling terms between elements of the heirarchy
% for r = 1:n_couplings
%     J = lower_indices(r) ;
%     K = upper_indices(r) ;
%     jk_coup = coupled_mode_indices(r) ;
%     j_coup = ceil(jk_coup/(M+1)) ; % get the bath index that is coupling J & K
%     k_coup = jk_coup - j_coup*(M+1) ; % Matsubara mode index, k=0 is non-Matsubara term
%     J_block = ((J-1)*(d_liou)+1):(J*d_liou) ;
%     K_block = ((K-1)*(d_liou)+1):(K*d_liou) ;
%     % add the coupling from the ADO deeper in the heirarchy 
%     % d/dt rho_J = ... -i [V_j , rho_K]
%     L_heom(J_block,K_block) = L_heom(J_block,K_block)-1.0i * V_comm{j_coup} ;
%     % add the coupling from the ADO shallower in the heirarchy
%     n_jk_coup = ado_indices(K,jk_coup) ;
%     c_jk_coup = cs(jk_coup) ;
%     L_heom(K_block,J_block) = L_heom(K_block,J_block) + (-1.0i*c_jk_coup*n_jk_coup) * V_L{j_coup} ...
%         + (1.0i*conj(c_jk_coup)*n_jk_coup) * V_R{j_coup} ; 
% end

% this uses the scaled couplings
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

% add a perturbative correction for the truncated ados
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
%     L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/sum(nus.*n_jks_term)) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;
% THIS IS THE main truncation used
    L_heom(J_block,J_block) = L_heom(J_block,J_block) -(1.0i/(ado_gammas(J)+nus(jk_term))) * V_comm{j_term}*((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
        + (1.0i*conj(cbar_jk_term)*(n_jk_term+1)) * V_R{j_term} ) ;

%     L_heom(J_block,J_block) = L_heom(J_block,J_block) -...
%         (1.0i) * V_comm{j_term}*((-L_sys-Xi+(ado_gammas(J)+nus(jk_term))*id_liou)\((-1.0i*c_jk_term*n_jk_term) * V_L{j_coup} ...
%         + (1.0i*conj(c_jk_term)*(n_jk_term+1)) * V_R{j_term} )) ;
end

end