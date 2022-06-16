function [K,scpt_junk] = constructHEOMSCPTCouplingOperator(full_system,heom_structure_blocks,heom_bath_info_blocks,block_coupling_info,heom_truncation_info)

% get the number of coupling terms between baths
n_couplings = size(full_system.block_coupling.coupled_blocks,1) ;

% first we need to generate the discretisations for each ET process
n_baths = []  ;
cs = {} ;
omegas = {} ;
weights = {} ;
lambda_hiT = {} ;
for r = 1:n_couplings
    n_baths(r) = numel(full_system.block_coupling.coupling_baths{r}) ;
    cs{r} = [] ;
    omegas{r} = [] ;
    weights{r} = [] ;
    lambda_hiT{r} = 0.0 ;
    for n = 1:n_baths(r)
        % construct the density of states for the spectral density
        bath_info = full_system.block_coupling.coupling_baths{r}{n} ;
        if bath_info.spectral_density == "debye"
            omega_D = bath_info.omega_D ;
            lambda_D = bath_info.lambda_D ;
            rho = @(x)(2.0/pi) * omega_D./(x.*x + omega_D*omega_D) ;
            % discretise the spectral density
            [omegas_bath,weights_bath] = discretiseSpecDenGL(bath_info.n_modes,rho) ;
            cs_bath = sqrt(2.0*lambda_D *weights_bath).*omegas_bath ;
            omegas{r}  = [omegas{r} , omegas_bath] ;
            cs{r} = [cs{r} , cs_bath] ;
            weights{r} = [weights{r},weights_bath * lambda_D] ;
        elseif bath_info.spectral_density == "BO"
            Omega = bath_info.Omega ;
            lambda = bath_info.lambda ;
            gamma = bath_info.gamma ;
            rho = @(x)(2.0/pi) * (gamma*Omega*Omega)./((x.*x-Omega*Omega).*(x.*x-Omega*Omega) + (gamma*gamma)*(x.*x)) ;
            % discretise the spectral density
            [omegas_bath,weights_bath] = discretiseSpecDenGL(bath_info.n_modes,rho) ;
            cs_bath = sqrt(2.0*lambda * weights_bath).*omegas_bath ;
            omegas{r}  = [omegas{r} , omegas_bath] ;
            cs{r} = [cs{r} , cs_bath] ;
            weights{r} = [weights{r},weights_bath * lambda] ;
        elseif (bath_info.spectral_density == "high temperature")
            lambda_hiT{r} = bath_info.lambda ;
        else
            % TODO: implement classical contribution to these correlation
            % functions.
            warning("This bath spectral density is not supported for the block coupling operator.")
        end
    end
end
t_maxs = full_system.block_coupling.t_maxs ;
n_ts = full_system.block_coupling.n_ts ;

% calculate the correlation function on the desired time grid
c_ts = cell([n_couplings,1]) ;
ts = cell([n_couplings,1]) ;
E_blocks = full_system.block_coupling.E_blocks ;
coupled_blocks = full_system.block_coupling.coupled_blocks ;
for r = 1:n_couplings
    t = linspace(0,t_maxs(r),n_ts(r)) ;
    n_A = coupled_blocks(r,1) ;
    n_B = coupled_blocks(r,2) ;
    Delta_E = E_blocks(n_A) - E_blocks(n_B) ;

    if (isempty(weights{r}))
        % <exp(-i H_B t) exp(+i H_A t)>_0
        c_A_ts = exp(1.0i * (Delta_E - lambda_hiT{r})*t - (lambda_hiT{r}/full_system.beta)*(t.*t)) ;
        % <exp(-i H_A t) exp(+i H_B t)>_B
        c_B_ts = exp(1.0i * (-Delta_E - lambda_hiT{r})*t - (lambda_hiT{r}/full_system.beta)*(t.*t)) ;
    else
        % <exp(-i H_B t) exp(+i H_A t)>_0
        c_A_ts = calculateCorrelationFunction(weights{r},omegas{r},t,full_system.beta,0)...
            .* exp(1.0i * (Delta_E - lambda_hiT{r})*t - (lambda_hiT{r}/full_system.beta)*(t.*t)) ;
        % <exp(-i H_A t) exp(+i H_B t)>_B
        c_B_ts = calculateCorrelationFunction(weights{r},omegas{r},t,full_system.beta,0)...
            .*exp(1.0i * (-Delta_E - lambda_hiT{r})*t - (lambda_hiT{r}/full_system.beta)*(t.*t)) ;
    end
    Gamma = full_system.block_coupling.coupling_matrices ;
    c_ts_r = [c_A_ts;c_B_ts] ;
    c_ts{r} = c_ts_r ;
    ts{r} = t ;
    int_c_t_A = trapz(t,c_A_ts) ;
    int_c_t_B = trapz(t,c_B_ts) ;
    fprintf('Rate term (forward):\n')
    Gamma{r}'*Gamma{r} * real(int_c_t_A) 
    fprintf('Shift term (init):\n')
    Gamma{r}'*Gamma{r} * imag(int_c_t_A) 
    fprintf('Rate term (back):\n')
    Gamma{r}*Gamma{r}' * real(int_c_t_B) 
    fprintf('Shift term (final):\n')
    Gamma{r}*Gamma{r}' * imag(int_c_t_B) 
    drawnow
end

% get general info for constructing the SCPT operator
n_blocks = numel(full_system.H_sys) ;
n_ados_blocks = [] ;
d_heoms = [] ;
d_lious = [] ;
d_hilbs = [] ;
block_indices = {};
for k = 1:n_blocks
    n_ados_blocks = [n_ados_blocks,heom_structure_blocks{k}.n_ados] ;
    d_heoms = [d_heoms,heom_structure_blocks{k}.d_heom] ;
    d_lious = [d_lious,heom_structure_blocks{k}.d_liou] ;
    d_hilbs = [d_hilbs,size(full_system.H_sys{k},1)] ;
    block_indices{k} = (1+sum(d_heoms(1:end-1))):sum(d_heoms(1:end)) ;
end
d = sum(d_heoms) ;
id_ados = speye(n_ados_blocks(1)) ;
fprintf('n_ado = %d\n',n_ados_blocks(1)) ;
K = sparse([],[],[],d,d) ;


% add in the coupling terms to the superoperator
if block_coupling_info.method == "truncated NZ"
    for r = 1:n_couplings
        n_A = coupled_blocks(r,1) ;
        n_B = coupled_blocks(r,2) ;
        % get the indices of the AB states included in the truncated NZ space
        if (numel(block_coupling_info.Gamma_cut_trunc) == 1)
            ado_inds = find(heom_structure_blocks{n_A}.ado_gammas < block_coupling_info.Gamma_cut_trunc) ;
        else 
            ado_inds = find(heom_structure_blocks{n_A}.ado_gammas < block_coupling_info.Gamma_cut_trunc(r)) ;
        end
        n_trunc_ado = length(ado_inds) ;
        fprintf('n_ado_trunc = %d.\n',n_trunc_ado) ;
        d_AB = d_hilbs(n_A) * d_hilbs(n_B) ;
        n_ados = n_ados_blocks(n_A) ;
        d_heom_AB = n_ados * d_AB ;
        d_hilb_A = d_hilbs(n_A) ;
        d_hilb_B = d_hilbs(n_B) ;
        d_heom_A = d_heoms(n_A) ;
        d_heom_B = d_heoms(n_B) ;
        id_hilb_A = speye(d_hilb_A) ;
        id_hilb_B = speye(d_hilb_B) ;
        inds_trunc = kron(ado_inds*d_AB,ones([d_AB,1])) - repmat(((d_AB-1):-1:0)',[n_trunc_ado,1]);
        inds_trunc_A = kron(ado_inds*d_hilb_A *d_hilb_A,ones([d_hilb_A *d_hilb_A,1])) - repmat(((d_hilb_A*d_hilb_A-1):-1:0)',[n_trunc_ado,1]);
        inds_trunc_B = kron(ado_inds*d_hilb_B *d_hilb_B,ones([d_hilb_B *d_hilb_B,1])) - repmat(((d_hilb_B*d_hilb_B-1):-1:0)',[n_trunc_ado,1]);
        d_trunc = length(inds_trunc) ;
        d_trunc_A = length(inds_trunc_A) ;
        d_trunc_B = length(inds_trunc_B) ;
        d_trunc_full = d_trunc_A + d_trunc_B ;
        % get the transformation matrix to truncate the space
        S_trunc = sparse([1:d_trunc]',[inds_trunc],ones([d_trunc,1]),d_trunc,d_heom_AB) ;
        S_trunc_A = sparse([1:d_trunc_A]',[inds_trunc_A],ones([d_trunc_A,1]),d_trunc_A,d_heom_A) ;
        S_trunc_B = sparse([1:d_trunc_B]',[inds_trunc_B],ones([d_trunc_B,1]),d_trunc_B,d_heom_B) ;
        S_trunc_full = sparse([],[],[],d_trunc_full, d) ;
        S_trunc_full(1:d_trunc_A,block_indices{n_A}) = S_trunc_A ;
        S_trunc_full(((d_trunc_A + 1):(d_trunc_A+d_trunc_B)),block_indices{n_B}) = S_trunc_B ;

        % construct the coherence blocks of the heom liouvillian
        H_sys_A = full_system.H_sys{n_A} ;
        H_sys_B = full_system.H_sys{n_B} ;
        V_As = heom_bath_info_blocks{n_A}.Vs ;
        V_Bs = heom_bath_info_blocks{n_B}.Vs ;
        heom_truncation_info_AB = heom_truncation_info ;
        if (isfield(block_coupling_info,'include_Xi_AB'))
            if (~block_coupling_info.include_Xi_AB)
                heom_truncation_info_AB.heom_termination = "none" ;
            end
        else 
            heom_truncation_info_AB.heom_termination = "none" ;
        end
        L_AB = constructHEOMABGenerator(H_sys_A,H_sys_B,V_As,V_Bs,heom_bath_info_blocks{n_A},heom_truncation_info_AB,heom_structure_blocks{n_A}) ;
        L_BA = constructHEOMABGenerator(H_sys_B,H_sys_A,V_Bs,V_As,heom_bath_info_blocks{n_A},heom_truncation_info_AB,heom_structure_blocks{n_A}) ;
        % truncate these
        L_AB_trunc = S_trunc *  L_AB * (S_trunc') ;
        L_BA_trunc = S_trunc *  L_BA * (S_trunc') ;

        % diagonalise the coherence blocks of the HEOM liouvillian
        [V_AB,Lambda_AB] = eig(full(L_AB_trunc),'vector') ;
        %     V_AB_inv = inv(V_AB) ;
        V_AB_inv = V_AB\eye(d_trunc) ;
        %     max(max(abs( V_AB*(Lambda_AB.*(V_AB_inv))-L_AB )))

        [V_BA,Lambda_BA] = eig(full(L_BA_trunc),'vector') ;
        %     V_AB_inv = inv(V_AB) ;
        V_BA_inv = V_BA\eye(d_trunc) ;
        %     max(max(abs( V_BA*(Lambda_BA.*(V_BA_inv))-L_BA )))

        % calculate some integrals
        t = ts{r} ; c_A_ts = c_ts{r}(1,:) ; c_B_ts = c_ts{r}(2,:) ;
        g_A_AB = trapz(t,conj(c_A_ts).*exp(Lambda_AB.*t),2) ;
        g_A_BA = trapz(t,c_A_ts.*exp(Lambda_BA.*t),2) ;
        g_B_AB = trapz(t,c_B_ts.*exp(Lambda_AB.*t),2) ;
        g_B_BA = trapz(t,conj(c_B_ts).*exp(Lambda_BA.*t),2) ;

        % construct some integrated effective propagators
        G_A_AB = V_AB * (g_A_AB.* V_AB_inv) ;
        G_A_BA = V_BA * (g_A_BA.* V_BA_inv) ;
        G_B_AB = V_AB * (g_B_AB.* V_AB_inv) ;
        G_B_BA = V_BA * (g_B_BA.* V_BA_inv) ;

        % construct the NZ2 rate operator in the truncated space
        id_ados_trunc = speye(n_trunc_ado) ;
        d_heom_A_trunc = d_hilb_A *  d_hilb_A * n_trunc_ado ;
        d_heom_B_trunc = d_hilb_B *  d_hilb_B * n_trunc_ado ;
        d_heom_trunc = d_heom_A_trunc + d_heom_B_trunc;
        K_trunc = zeros([d_heom_trunc,d_heom_trunc]) ;

        % construct the AA part of the NZ2 rate operator
        Gamma = full_system.block_coupling.coupling_matrices{r} ;
        Gamma_L_AA = kron(id_ados_trunc,kron(Gamma',id_hilb_A)) ;
        Gamma_L_BA = kron(id_ados_trunc,kron(Gamma,id_hilb_A)) ;
        Gamma_R_AA = kron(id_ados_trunc,kron(id_hilb_A,transpose(Gamma))) ;
        Gamma_R_AB = kron(id_ados_trunc,kron(id_hilb_A,transpose(Gamma'))) ;

        K_trunc(1:d_heom_A_trunc,1:d_heom_A_trunc) = -Gamma_L_BA * G_A_BA  * Gamma_L_AA  - Gamma_R_AB * G_A_AB * Gamma_R_AA ;

        % construct the BB part
        Gamma_L_BB = kron(id_ados_trunc,kron(Gamma,id_hilb_B)) ;
        Gamma_L_AB = kron(id_ados_trunc,kron(Gamma',id_hilb_B)) ;
        Gamma_R_BB = kron(id_ados_trunc,kron(id_hilb_B,transpose(Gamma'))) ;
        Gamma_R_BA = kron(id_ados_trunc,kron(id_hilb_B,transpose(Gamma))) ;
        K_trunc((d_heom_A_trunc+1):(d_heom_A_trunc+d_heom_B_trunc),(d_heom_A_trunc+1):(d_heom_A_trunc+d_heom_B_trunc)) = ...
            -Gamma_L_AB * G_B_AB * Gamma_L_BB  - Gamma_R_BA * G_B_BA * Gamma_R_BB;

        % construct BA part
        K_trunc((d_heom_A_trunc+1):(d_heom_A_trunc+d_heom_B_trunc),1:d_heom_A_trunc) = ...
            Gamma_R_BA * G_A_BA  * Gamma_L_AA  + Gamma_L_AB * G_A_AB * Gamma_R_AA ;

        % construct the AB part
        K_trunc(1:d_heom_A_trunc,(d_heom_A_trunc+1):(d_heom_A_trunc+d_heom_B_trunc)) = ...
            Gamma_R_AB * G_B_AB * Gamma_L_BB  + Gamma_L_BA * G_B_BA * Gamma_R_BB ;

        % construct K for the remaining ados
        % construct the L_sys operators acting on AB and BA spaces
        L_sys_AB = -1.0i * kron(H_sys_A,id_hilb_B) +1.0i *kron(id_hilb_A,transpose(H_sys_B)) ;
        L_sys_BA = -1.0i * kron(H_sys_B,id_hilb_A) +1.0i *kron(id_hilb_B,transpose(H_sys_A)) ;
        [P_AB,Lambda_AB] = eig(full(L_sys_AB),'vector') ;
        P_AB_inv = P_AB\eye(d_hilb_A*d_hilb_B) ;
        [P_BA,Lambda_BA] = eig(full(L_sys_BA),'vector') ;
        P_BA_inv = P_BA\eye(d_hilb_A*d_hilb_B) ;
        Pi_AB = kron(id_ados,P_AB) ;
        Pi_AB_inv = kron(id_ados,P_AB_inv) ;
        Pi_BA = kron(id_ados,P_BA) ;
        Pi_BA_inv = kron(id_ados,P_BA_inv) ;
        % construct some fourier transforms
        %         n_rep = n_ados ;
        %         g_A_AB_0 = repmat(trapz(t,conj(c_A_ts).*exp(Lambda_AB.*t),2),[n_rep,1]) ;
        %         g_A_BA_0 = repmat(trapz(t,c_A_ts.*exp(Lambda_BA.*t),2),[n_rep,1]) ;
        %         g_B_AB_0 = repmat(trapz(t,c_B_ts.*exp(Lambda_AB.*t),2),[n_rep,1]) ;
        %         g_B_BA_0 = repmat(trapz(t,conj(c_B_ts).*exp(Lambda_BA.*t),2),[n_rep,1]) ;
        %         Lambda_AB = repmat(Lambda_AB,[n_rep,1]) ;
        %         Lambda_BA = repmat(Lambda_BA,[n_rep,1]) ;

        n_rep = 1 ;
        Lambda_AB = kron(ones([n_ados,1]),Lambda_AB) - kron((heom_structure_blocks{n_A}.ado_gammas),ones([d_AB,1]));
        Lambda_BA = kron(ones([n_ados,1]),Lambda_BA) - kron((heom_structure_blocks{n_A}.ado_gammas),ones([d_AB,1]));
        g_A_AB_0 = repmat(trapz(t,conj(c_A_ts).*exp(Lambda_AB.*t),2),[n_rep,1]) ;
        g_A_BA_0 = repmat(trapz(t,c_A_ts.*exp(Lambda_BA.*t),2),[n_rep,1]) ;
        g_B_AB_0 = repmat(trapz(t,c_B_ts.*exp(Lambda_AB.*t),2),[n_rep,1]) ;
        g_B_BA_0 = repmat(trapz(t,conj(c_B_ts).*exp(Lambda_BA.*t),2),[n_rep,1]) ;

        % find the indices of the AB ados in the complementary space
        inds_comp = 1:d_heom_AB ;
        is_trunc = ismember(inds_comp,inds_trunc) ;
        inds_comp = inds_comp(~is_trunc) ;
        g_A_BA_0(is_trunc) = 0 ;
        g_A_BA_0(is_trunc) = 0 ;
        g_B_AB_0(is_trunc) = 0 ;
        g_B_BA_0(is_trunc) = 0 ;
        % construct the remaining portion on the operator
        G_A_BA_0 = Pi_BA * ((repmat(g_A_BA_0,[1,1]).* Pi_BA_inv)) ;
        G_A_AB_0 = Pi_AB * ((repmat(g_A_AB_0,[1,1]).* Pi_AB_inv)) ;
        G_B_BA_0 = Pi_BA * ((repmat(g_B_BA_0,[1,1]).* Pi_BA_inv)) ;
        G_B_AB_0 = Pi_AB * ((repmat(g_B_AB_0,[1,1]).* Pi_AB_inv)) ;

        % construct Gamma operators
        Gamma = full_system.block_coupling.coupling_matrices{r} ;
        if (numel(block_coupling_info.Gamma_cut_trunc) == 1)
            proj_comp_ados = spdiags(heom_structure_blocks{n_A}.ado_gammas >= block_coupling_info.Gamma_cut_trunc,[0],n_ados,n_ados);
        else 
            proj_comp_ados = spdiags(heom_structure_blocks{n_A}.ado_gammas >= block_coupling_info.Gamma_cut_trunc(r),[0],n_ados,n_ados);
        end
        
        Gamma_L_AA = kron(proj_comp_ados,kron(Gamma',id_hilb_A)) ;
        Gamma_L_BA = kron(proj_comp_ados,kron(Gamma,id_hilb_A)) ;
        Gamma_R_AA = kron(proj_comp_ados,kron(id_hilb_A,transpose(Gamma))) ;
        Gamma_R_AB = kron(proj_comp_ados,kron(id_hilb_A,transpose(Gamma'))) ;
        Gamma_L_BB = kron(proj_comp_ados,kron(Gamma,id_hilb_B)) ;
        Gamma_L_AB = kron(proj_comp_ados,kron(Gamma',id_hilb_B)) ;
        Gamma_R_BB = kron(proj_comp_ados,kron(id_hilb_B,transpose(Gamma'))) ;
        Gamma_R_BA = kron(proj_comp_ados,kron(id_hilb_B,transpose(Gamma))) ;


        % BB part
        K(block_indices{n_B},block_indices{n_B}) = K(block_indices{n_B},block_indices{n_B})  ...
            - Gamma_L_AB * (G_B_AB_0) * Gamma_L_BB  - Gamma_R_BA * (G_B_BA_0) * Gamma_R_BB;
        % construct the AB part
        K(block_indices{n_A},block_indices{n_B}) = K(block_indices{n_A},block_indices{n_B}) ...
            + Gamma_R_AB * (G_B_AB_0) * Gamma_L_BB  + Gamma_L_BA * (G_B_BA_0) * Gamma_R_BB ;
        % AA part
        K(block_indices{n_A},block_indices{n_A}) = K(block_indices{n_A},block_indices{n_A}) ...
            - Gamma_L_BA * (G_A_BA_0)  * Gamma_L_AA  - Gamma_R_AB * (G_A_AB_0) * Gamma_R_AA ;

        % construct BA part
        K(block_indices{n_B},block_indices{n_A}) = K(block_indices{n_B},block_indices{n_A}) ...
            + Gamma_R_BA * (G_A_BA_0)  * Gamma_L_AA  + Gamma_L_AB * (G_A_AB_0) * Gamma_R_AA ;


        if (n_trunc_ado == n_ados)
            K = K + S_trunc_full' * full(K_trunc) * S_trunc_full ;
        else
            K = K + S_trunc_full' * sparse(K_trunc) * S_trunc_full ;
        end
    end

end
drawnow('update')


% add in the approximate radiative decay part if specified
if (isfield(full_system.block_coupling,'coupled_blocks_radiative'))
    fprintf('Adding radiative coupling terms.\n')
    n_couplings = size(full_system.block_coupling.coupled_blocks_radiative,1) ;
    for r = 1:n_couplings

        n_A = full_system.block_coupling.coupled_blocks_radiative(r,1) ;
        n_B = full_system.block_coupling.coupled_blocks_radiative(r,2) ;
        d_AB = d_hilbs(n_A) * d_hilbs(n_B) ;
        n_ados = n_ados_blocks(n_A) ;
        d_heom_AB = n_ados * d_AB ;
        d_hilb_A = d_hilbs(n_A) ;
        d_hilb_B = d_hilbs(n_B) ;
        d_heom_A = d_heoms(n_A) ;
        d_heom_B = d_heoms(n_B) ;
        id_hilb_A = speye(d_hilb_A) ;
        id_hilb_B = speye(d_hilb_B) ;
        % construct K for the remaining ados
        % construct the L_sys operators acting on AB and BA spaces
        H_sys_A = full_system.H_sys{n_A} ;
        H_sys_B = full_system.H_sys{n_B} ;
        L_sys_AB = -1.0i * kron(H_sys_A,id_hilb_B) +1.0i *kron(id_hilb_A,transpose(H_sys_B)) ;
        L_sys_BA = -1.0i * kron(H_sys_B,id_hilb_A) +1.0i *kron(id_hilb_B,transpose(H_sys_A)) ;
        [P_AB,Lambda_AB] = eig(full(L_sys_AB),'vector') ;
        Omega_AB = imag(Lambda_AB) ;
        P_AB_inv = P_AB\eye(d_hilb_A*d_hilb_B) ;
        [P_BA,Lambda_BA] = eig(full(L_sys_BA),'vector') ;
        Omega_BA = imag(Lambda_BA ) ;
        P_BA_inv = P_BA\eye(d_hilb_A*d_hilb_B) ;
        Pi_AB = kron(id_ados,P_AB) ;
        Pi_AB_inv = kron(id_ados,P_AB_inv) ;
        Pi_BA = kron(id_ados,P_BA) ;
        Pi_BA_inv = kron(id_ados,P_BA_inv) ;
        % construct some fourier transforms - check signs here
        Delta_E_AB = E_blocks(n_A) - E_blocks(n_B) ;
        Delta_E_BA = E_blocks(n_B) - E_blocks(n_A) ;
        g_B_BA_0 = pi * (-(Omega_BA-Delta_E_BA).^3) ;
        g_B_BA_0(g_B_BA_0<0) = 0 ;
        g_B_BA_0 = repmat(g_B_BA_0,[n_ados,1]) ;
        g_B_AB_0 = g_B_BA_0 ;

        g_A_BA_0 = pi * ((Omega_BA-Delta_E_BA).^3) ;
        g_A_BA_0(g_A_BA_0<0) = 0 ;
        g_A_BA_0 = repmat(g_A_BA_0,[n_ados,1]) ;
        g_A_AB_0 = g_A_BA_0 ;


        % construct the remaining portion on the operator
        G_A_BA_0 = Pi_BA * ((repmat(g_A_BA_0,[1,1]).* Pi_BA_inv)) ;
        G_A_AB_0 = Pi_AB * ((repmat(g_A_AB_0,[1,1]).* Pi_AB_inv)) ;
        G_B_BA_0 = Pi_BA * ((repmat(g_B_BA_0,[1,1]).* Pi_BA_inv)) ;
        G_B_AB_0 = Pi_AB * ((repmat(g_B_AB_0,[1,1]).* Pi_AB_inv)) ;

        n_trans_ops = numel(full_system.block_coupling.coupling_matrices_radiative{r}) ;
        for alpha = 1:n_trans_ops
            % construct Gamma operators
            Gamma = full_system.block_coupling.coupling_matrices_radiative{r}{alpha} ;
            proj_comp_ados = speye(n_ados) ;
            Gamma_L_AA = kron(proj_comp_ados,kron(Gamma',id_hilb_A)) ;
            Gamma_L_BA = kron(proj_comp_ados,kron(Gamma,id_hilb_A)) ;
            Gamma_R_AA = kron(proj_comp_ados,kron(id_hilb_A,transpose(Gamma))) ;
            Gamma_R_AB = kron(proj_comp_ados,kron(id_hilb_A,transpose(Gamma'))) ;
            Gamma_L_BB = kron(proj_comp_ados,kron(Gamma,id_hilb_B)) ;
            Gamma_L_AB = kron(proj_comp_ados,kron(Gamma',id_hilb_B)) ;
            Gamma_R_BB = kron(proj_comp_ados,kron(id_hilb_B,transpose(Gamma'))) ;
            Gamma_R_BA = kron(proj_comp_ados,kron(id_hilb_B,transpose(Gamma))) ;

            % BB part
            K(block_indices{n_B},block_indices{n_B}) = K(block_indices{n_B},block_indices{n_B})  ...
                - Gamma_L_AB * (G_B_AB_0) * Gamma_L_BB  - Gamma_R_BA * (G_B_BA_0) * Gamma_R_BB;
            % construct the AB part
            K(block_indices{n_A},block_indices{n_B}) = K(block_indices{n_A},block_indices{n_B}) ...
                + Gamma_R_AB * (G_B_AB_0) * Gamma_L_BB  + Gamma_L_BA * (G_B_BA_0) * Gamma_R_BB ;
            % AA part
            K(block_indices{n_A},block_indices{n_A}) = K(block_indices{n_A},block_indices{n_A}) ...
                - Gamma_L_BA * (G_A_BA_0)  * Gamma_L_AA  - Gamma_R_AB * (G_A_AB_0) * Gamma_R_AA ;

            % construct BA part
            K(block_indices{n_B},block_indices{n_A}) = K(block_indices{n_B},block_indices{n_A}) ...
                + Gamma_R_BA * (G_A_BA_0)  * Gamma_L_AA  + Gamma_L_AB * (G_A_AB_0) * Gamma_R_AA ;
        end
    end
end

scpt_junk = {ts,c_ts} ;


end