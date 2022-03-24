function [K,c_ts,ts] = constructTC2ABOperator(AB_coupling_info, H_sys_A, H_sys_B,d_heom_A,d_heom_B,heom_structure_A,heom_structure_B,V_As,V_Bs,heom_bath_info,heom_truncation_info)

% first we need to generate the correlation functions for each bath
n_baths = numel(AB_coupling_info.baths) ;
cs = [] ;
omegas = [] ;
weights = [] ;
for n = 1:n_baths
    % construct the density of states for the spectral density
    if AB_coupling_info.baths{n}.spectral_density == "debye"
        omega_D = AB_coupling_info.baths{n}.omega_D ;
        lambda_D = AB_coupling_info.baths{n}.lambda_D ;
        rho = @(x)(2.0/pi) * omega_D./(x.*x + omega_D*omega_D) ;
        % discretise the spectral density
        [omegas_bath,weights_bath] = discretiseSpecDenGL(AB_coupling_info.n_modes,rho) ;
        cs_bath = sqrt(2.0*lambda_D *weights_bath).*omegas_bath ;
        omegas  = [omegas , omegas_bath] ;
        cs = [cs , cs_bath] ;
        weights = [weights,weights_bath * lambda_D] ;
    end

end

% calculate the correlation function on the desired time grid
ts = linspace(0,AB_coupling_info.t_max,AB_coupling_info.n_t) ;
% <exp(-i H_B t) exp(+i H_A t)>_A
c_A_ts = calculateCorrelationFunction(weights,omegas,ts,AB_coupling_info.beta,AB_coupling_info.Delta_E_AB) ;
% <exp(-i H_A t) exp(+i H_B t)>_B
c_B_ts = calculateCorrelationFunction(weights,omegas,ts,AB_coupling_info.beta,-AB_coupling_info.Delta_E_AB) ;
c_ts = [c_A_ts;c_B_ts] ;
% get dimensions of the system spaces
d_hilb_A = size(H_sys_A,1) ;
d_hilb_B = size(H_sys_B,1) ;
d = d_heom_A + d_heom_B ;
d_liou_A = d_hilb_A^2 ;
d_liou_B = d_hilb_B^2 ;
K = sparse([],[],[],d,d);
id_hilb_A = speye(d_hilb_A) ;
id_hilb_B = speye(d_hilb_B) ;
n_ados = d_heom_A / d_liou_A ;
id_ados = speye(n_ados) ;

% the simplified version simply sets
if AB_coupling_info.method == "simplified"
    int_c_A = trapz(ts,c_A_ts) ;
    int_c_B = trapz(ts,c_B_ts) ;
    Gamma = AB_coupling_info.coupling_matrix ;
    % construct the AA transfer operator
    T_AA = -int_c_A * Gamma*(Gamma') ;
    % construct the AA rate superoperator K_AA sigma = T_AA sigma + sigma T_AA^dag
    K_AA_liou = kron(T_AA,id_hilb_A) + kron(id_hilb_A,conj(T_AA)) ;
    K(1:d_heom_A,1:d_heom_A) = kron(id_ados,K_AA_liou) ;

    % construct T_BB
    T_BB = -int_c_B * ((Gamma')*(Gamma)) ;
    K_BB_liou = kron(T_BB,id_hilb_B) + kron(id_hilb_B,conj(T_BB)) ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_BB_liou) ;

    % construct the B<-A term
    K_BA_liou = (int_c_A + conj(int_c_A))*kron(Gamma',transpose(Gamma)) ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = kron(id_ados,K_BA_liou) ;
    % construct the A<-B term
    K_AB_liou = (int_c_B + conj(int_c_B))*kron(Gamma,conj(Gamma)) ;
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_AB_liou) ;
elseif AB_coupling_info.method == "include H_sys"
    % calculate the eigenstates of the system Hamiltonians
    [psi_As,E_As] = eig(H_sys_A,'vector') ;
    [psi_Bs,E_Bs] = eig(H_sys_B,'vector') ;


    % calculate the coupling matrix in the system energy eigenbasis
    Gamma_E = psi_As' * AB_coupling_info.coupling_matrix * psi_Bs ;
    Delta_E_AB_sys = E_Bs - E_As' ;

    G_A = zeros(d_hilb_B,d_hilb_A) ;
    G_B = zeros(d_hilb_A,d_hilb_B) ;
    for n_A = 1:d_hilb_A
        for n_B = 1:d_hilb_B
            G_A(n_B,n_A) = trapz(ts,c_A_ts.*exp(-1.0i*(E_Bs(n_B)-E_As(n_A))*ts)) ;
            % need to check this G_B result
            G_B(n_A,n_B) = trapz(ts,(c_B_ts).*exp(+1.0i*(E_Bs(n_B)-E_As(n_A))*ts)) ;
        end
    end

    % construct the AA rate superoperator K_AA sigma = T_AA sigma + sigma T_AA^dag
    T_AA_E = -Gamma_E * (G_A .* Gamma_E') ;
    T_AA = psi_As*T_AA_E*(psi_As')  ;
    K_AA_liou = kron(T_AA,id_hilb_A) + kron(id_hilb_A,conj(T_AA)) ;
    K(1:d_heom_A,1:d_heom_A) = kron(id_ados,K_AA_liou) ;
    %     % construct T_BB
    T_BB_E =  -Gamma_E' * (G_B.*Gamma_E) ;
    T_BB = psi_Bs*T_BB_E*(psi_Bs')  ;
    K_BB_liou = kron(T_BB,id_hilb_B) + kron(id_hilb_B,conj(T_BB)) ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_BB_liou) ;
    % construct the A->B term
    K_BA_liou = kron(psi_Bs*(G_A.*Gamma_E')*(psi_As'),transpose(psi_As*Gamma_E*(psi_Bs')))...
        + kron(psi_Bs*Gamma_E'*(psi_As'),conj( (psi_Bs*(G_A.*Gamma_E')*(psi_As')) ) );
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = kron(id_ados,K_BA_liou) ;
    %     % construct the A<-B term
    K_AB_liou = kron(psi_As*(G_B.*Gamma_E)*(psi_Bs'),transpose(psi_Bs*Gamma_E'*(psi_As')))...
        + kron(psi_As*Gamma_E*(psi_Bs'),conj( (psi_As*(G_B.*Gamma_E)*(psi_Bs')) ) );
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_AB_liou) ;

elseif AB_coupling_info.method == "include H_sys NZ"
    % calculate the eigenstates of the system Hamiltonians
    [psi_As,E_As] = eig(H_sys_A,'vector') ;
    [psi_Bs,E_Bs] = eig(H_sys_B,'vector') ;


    % calculate the coupling matrix in the system energy eigenbasis
    Gamma_E = psi_As' * AB_coupling_info.coupling_matrix * psi_Bs ;
    Gamma = AB_coupling_info.coupling_matrix ;
    Delta_E_AB_sys = E_Bs - E_As' ;
    int_c_A = trapz(ts,c_A_ts) ;
    int_c_B = trapz(ts,c_B_ts) ;
    G_A = zeros(d_hilb_B,d_hilb_A) ;
    G_B = zeros(d_hilb_A,d_hilb_B) ;
    for n_A = 1:d_hilb_A
        for n_B = 1:d_hilb_B
            G_A(n_B,n_A) = trapz(ts,c_A_ts.*exp(-1.0i*(E_Bs(n_B)-E_As(n_A))*ts)) ;
            % need to check this G_B result
            G_B(n_A,n_B) = trapz(ts,(c_B_ts).*exp(+1.0i*(E_Bs(n_B)-E_As(n_A))*ts)) ;
        end
    end
    % site to energy eigenbasis transform
    U_A = kron(psi_As,conj(psi_As))' ;
    U_B = kron(psi_Bs,conj(psi_Bs))' ;

    % K_AA term
    K_AA_liou = zeros([d_liou_A,d_liou_A]) ;
    for n=1:d_hilb_A
        P_An = zeros([d_hilb_A,d_hilb_A]) ;
        P_An(n,n) = 1 ;
        K_AA_liou = K_AA_liou - kron(Gamma_E * ( G_A(:,n).* Gamma_E' ),transpose(P_An))...
            - kron(P_An',conj(Gamma_E * (G_A(:,n) .* Gamma_E' )));
    end
    K_AA_liou = U_A' * K_AA_liou * U_A ;


    K(1:d_heom_A,1:d_heom_A) = kron(id_ados,K_AA_liou) ;
    % K_BB term
    K_BB_liou = zeros([d_liou_B,d_liou_B]) ;
    for n=1:d_hilb_B
        P_Bn = zeros([d_hilb_B,d_hilb_B]) ;
        P_Bn(n,n) = 1 ;
        K_BB_liou = K_BB_liou - kron(Gamma_E' * (Gamma_E .* G_B(:,n)),transpose(P_Bn))...
            - kron(P_Bn',conj(Gamma_E' * (Gamma_E .* G_B(:,n))));
    end
    K_BB_liou = U_B' * K_BB_liou * U_B ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_BB_liou) ;
    % K_BA term
    K_BA_liou = zeros([d_liou_B,d_liou_A]) ;
    for n=1:d_hilb_A
        Gamma_E_n = zeros([d_hilb_A,d_hilb_B]) ;
        Gamma_E_n(n,:) = Gamma_E(n,:) ;
        K_BA_liou = K_BA_liou + kron((Gamma_E' .* G_A(:,n)),transpose(Gamma_E_n))...
            + kron(Gamma_E_n',conj(Gamma_E' .* G_A(:,n)));
    end
    K_BA_liou = U_B' * K_BA_liou * U_A ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = kron(id_ados,K_BA_liou) ;
    % K_AB term
    K_AB_liou = zeros([d_liou_A,d_liou_B]) ;
    for n=1:d_hilb_B
        Gamma_E_n = zeros([d_hilb_A,d_hilb_B]) ;
        Gamma_E_n(:,n) = Gamma_E(:,n) ;
        K_AB_liou = K_AB_liou + kron((Gamma_E .* G_B(:,n)),transpose(Gamma_E_n'))...
            + kron(Gamma_E_n,conj(Gamma_E .* G_B(:,n)));
    end
    % construct the A<-B term

    K_AB_liou = U_A' * K_AB_liou * U_B ;
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = kron(id_ados,K_AB_liou) ;

elseif AB_coupling_info.method == "full NZ"
    % construct the coherence blocks of the heom liouvillian
    L_AB = constructHEOMABGenerator(H_sys_A,H_sys_B,V_As,V_Bs,heom_bath_info,heom_truncation_info,heom_structure_A) ;
    L_BA = constructHEOMABGenerator(H_sys_B,H_sys_A,V_Bs,V_As,heom_bath_info,heom_truncation_info,heom_structure_A) ;

    % diagonalise the coherence blocks of the HEOM liouvillian
    [V_AB,Lambda_AB] = eig(full(L_AB),'vector') ;
    %     V_AB_inv = inv(V_AB) ;
    V_AB_inv = V_AB\eye(d_hilb_A*d_hilb_B*n_ados) ;
    %     max(max(abs( V_AB*(Lambda_AB.*(V_AB_inv))-L_AB )))

    [V_BA,Lambda_BA] = eig(full(L_BA),'vector') ;
    %     V_AB_inv = inv(V_AB) ;
    V_BA_inv = V_BA\eye(d_hilb_A*d_hilb_B*n_ados) ;
    %     max(max(abs( V_BA*(Lambda_BA.*(V_BA_inv))-L_BA )))

    % calculate some integrals
    g_A_AB = trapz(ts,conj(c_A_ts).*exp(Lambda_AB.*ts),2) ;
    g_A_BA = trapz(ts,c_A_ts.*exp(Lambda_BA.*ts),2) ;
    g_B_AB = trapz(ts,c_B_ts.*exp(Lambda_AB.*ts),2) ;
    g_B_BA = trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2) ;
    %     g_B_AB = trapz(ts,(c_B_ts).*exp(Lambda_AB.*ts),2) ;
    %     g_B_BA = trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2) ;

    % construct some integrated effective propagators
    G_A_AB = V_AB * (g_A_AB.* V_AB_inv) ;
    G_A_BA = V_BA * (g_A_BA.* V_BA_inv) ;
    G_B_AB = V_AB * (g_B_AB.* V_AB_inv) ;
    G_B_BA = V_BA * (g_B_BA.* V_BA_inv) ;

    % construct the NZ2 rate operator
    K = zeros([d,d]) ;

    % construct the AA part of the NZ2 rate operator
    Gamma = AB_coupling_info.coupling_matrix ;
    Gamma_L_AA = kron(id_ados,kron(Gamma',id_hilb_A)) ;
    Gamma_L_BA = kron(id_ados,kron(Gamma,id_hilb_A)) ;
    Gamma_R_AA = kron(id_ados,kron(id_hilb_A,transpose(Gamma))) ;
    Gamma_R_AB = kron(id_ados,kron(id_hilb_A,transpose(Gamma'))) ;

    K(1:d_heom_A,1:d_heom_A) = -Gamma_L_BA * G_A_BA  * Gamma_L_AA  - Gamma_R_AB * G_A_AB * Gamma_R_AA ;

    % construct the BB part
    Gamma_L_BB = kron(id_ados,kron(Gamma,id_hilb_B)) ;
    Gamma_L_AB = kron(id_ados,kron(Gamma',id_hilb_B)) ;
    Gamma_R_BB = kron(id_ados,kron(id_hilb_B,transpose(Gamma'))) ;
    Gamma_R_BA = kron(id_ados,kron(id_hilb_B,transpose(Gamma))) ;
    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        -Gamma_L_AB * G_B_AB * Gamma_L_BB  - Gamma_R_BA * G_B_BA * Gamma_R_BB;

    % construct BA part
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = ...
        Gamma_R_BA * G_A_BA  * Gamma_L_AA  + Gamma_L_AB * G_A_AB * Gamma_R_AA ;

    % construct the AB part
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        Gamma_R_AB * G_B_AB * Gamma_L_BB  + Gamma_L_BA * G_B_BA * Gamma_R_BB ;

elseif (AB_coupling_info.method == "first-order phonon NZ" || AB_coupling_info.method == "first-order phonon NZ 2")
    % construct the coherence blocks of the heom liouvillian
    L_phonon_AB = constructHEOMABGenerator(0*H_sys_A,0*H_sys_B,V_As,V_Bs,heom_bath_info,heom_truncation_info,heom_structure_A) ;
    L_phonon_BA = constructHEOMABGenerator(0*H_sys_B,0*H_sys_A,V_Bs,V_As,heom_bath_info,heom_truncation_info,heom_structure_A) ;
    if (AB_coupling_info.method == "first-order phonon NZ 2")
        L_phonon_AB = L_phonon_AB + kron(spdiags([heom_structure_A.ado_gammas],[0],n_ados,n_ados),speye(d_hilb_A*d_hilb_B)) ;
        L_phonon_BA = L_phonon_BA + kron(spdiags([heom_structure_A.ado_gammas],[0],n_ados,n_ados),speye(d_hilb_A*d_hilb_B)) ;

    end

    % construct Gamma operators
    Gamma = AB_coupling_info.coupling_matrix ;
    Gamma_L_AA = kron(id_ados,kron(Gamma',id_hilb_A)) ;
    Gamma_L_BA = kron(id_ados,kron(Gamma,id_hilb_A)) ;
    Gamma_R_AA = kron(id_ados,kron(id_hilb_A,transpose(Gamma))) ;
    Gamma_R_AB = kron(id_ados,kron(id_hilb_A,transpose(Gamma'))) ;
    Gamma_L_BB = kron(id_ados,kron(Gamma,id_hilb_B)) ;
    Gamma_L_AB = kron(id_ados,kron(Gamma',id_hilb_B)) ;
    Gamma_R_BB = kron(id_ados,kron(id_hilb_B,transpose(Gamma'))) ;
    Gamma_R_BA = kron(id_ados,kron(id_hilb_B,transpose(Gamma))) ;

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
    L_tilde_phonon_AB = Pi_AB_inv * L_phonon_AB * Pi_AB ;
    L_tilde_phonon_BA = Pi_BA_inv * L_phonon_BA * Pi_BA ;
    d_AB = d_hilb_A * d_hilb_B ;
    % construct some fourier transforms
    if (AB_coupling_info.method == "first-order phonon NZ")
        g_A_AB = repmat(trapz(ts,conj(c_A_ts).*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_A_BA = repmat(trapz(ts,c_A_ts.*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_B_AB = repmat(trapz(ts,c_B_ts.*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_B_BA = repmat(trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_A_AB_1 = repmat(trapz(ts,ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_A_BA_1 = repmat(trapz(ts,ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_B_AB_1 = repmat(trapz(ts,ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_B_BA_1 = repmat(trapz(ts,ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        Lambda_AB = repmat(Lambda_AB,[n_ados,1]) ;
        Lambda_BA = repmat(Lambda_BA,[n_ados,1]) ;
    elseif (AB_coupling_info.method == "first-order phonon NZ 2")
        Lambda_AB = repmat(Lambda_AB,[n_ados,1]) - kron(heom_structure_A.ado_gammas,ones([d_hilb_A *d_hilb_B,1]));
        Lambda_BA = repmat(Lambda_BA,[n_ados,1]) - kron(heom_structure_A.ado_gammas,ones([d_hilb_A *d_hilb_B,1]));
        g_A_AB = repmat(trapz(ts,conj(c_A_ts).*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_A_BA = repmat(trapz(ts,c_A_ts.*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_B_AB = repmat(trapz(ts,c_B_ts.*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_B_BA = repmat(trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_A_AB_1 = repmat(trapz(ts,ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_A_BA_1 = repmat(trapz(ts,ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_B_AB_1 = repmat(trapz(ts,ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_B_BA_1 = repmat(trapz(ts,ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
    end

    % K_AA term
    G_A_BA_0 = Pi_BA * (repmat(g_A_BA,[1,1]).* Pi_BA_inv) ;
    G_A_AB_0 = Pi_AB * (repmat(g_A_AB,[1,1]).* Pi_AB_inv) ;
    %     j_A_BA = ((-repmat(g_A_BA,[1,1]) + transpose(repmat(g_A_BA,[1,1])))./(-repmat(Lambda_BA,[1,1])+transpose(repmat(Lambda_BA,[1,1])))) ;
    %     [rows,cols] = find(isnan(j_A_BA)) ;
    %     for r = 1:numel(rows)
    %         j_A_BA(rows(r),cols(r)) = g_A_BA_1(rows(r)) ;
    %     end
    %     %     j_A_BA(isnan(j_A_BA)) = g_A_BA_1(isnan(j_A_BA)) ;
    %     j_A_AB = ((-repmat(g_A_AB,[1,1]) + transpose(repmat(g_A_AB,[1,1])))./(-repmat(Lambda_AB,[1,1])+transpose(repmat(Lambda_AB,[1,1])))) ;
    %     [rows,cols] = find(isnan(j_A_AB)) ;
    %     for r = 1:numel(rows)
    %         j_A_AB(rows(r),cols(r)) = g_A_AB_1(rows(r)) ;
    %     end
    %     j_A_AB(isnan(j_A_AB)) = g_A_AB_1(isnan(j_A_AB));
    %     G_A_BA_1 = Pi_BA*((Pi_BA_inv * L_phonon_BA * Pi_BA).* ((repmat(g_A_BA,[n_ados,1]) - transpose(repmat(g_A_BA,[n_ados,1])))./(-repmat(Lambda_BA,[n_ados,1])+transpose(repmat(Lambda_BA,[n_ados,1])))))*Pi_BA_inv ;
    %     G_A_AB_1 = Pi_AB*((Pi_AB_inv * L_phonon_AB * Pi_AB).* ((repmat(g_A_AB,[n_ados,1]) - transpose(repmat(g_A_AB,[n_ados,1])))./(-repmat(Lambda_AB,[n_ados,1])+transpose(repmat(Lambda_AB,[n_ados,1])))))*Pi_AB_inv ;
    %     G_A_BA_1 = L_phonon_BA *g_A_BA_1 ;
    %     G_A_AB_1 = L_phonon_AB *g_A_AB_1 ;
    %     G_A_BA_1 = Pi_BA*((Pi_BA_inv * L_phonon_BA * Pi_BA).* j_A_BA)*Pi_BA_inv ;
    %     G_A_AB_1 = Pi_AB*((Pi_AB_inv * L_phonon_AB * Pi_AB).* j_A_AB)*Pi_AB_inv ;


    G_A_BA_1 = L_tilde_phonon_BA ;
    [rows,cols] = find(G_A_BA_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            G_A_BA_1(rows(r),cols(r)) = G_A_BA_1(rows(r),cols(r)) * g_A_BA_1(n) ;
        else
            G_A_BA_1(rows(r),cols(r)) = G_A_BA_1(rows(r),cols(r)) * (g_A_BA(n)-g_A_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
        end
    end
    G_A_BA_1 = Pi_BA * G_A_BA_1 * Pi_BA_inv ;

    G_A_AB_1 = L_tilde_phonon_AB ;
    [rows,cols] = find(G_A_AB_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            G_A_AB_1(rows(r),cols(r)) = G_A_AB_1(rows(r),cols(r)) * g_A_AB_1(n) ;
        else
            G_A_AB_1(rows(r),cols(r)) = G_A_AB_1(rows(r),cols(r)) * (g_A_AB(n)-g_A_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
        end
    end
    G_A_AB_1 = Pi_AB * G_A_AB_1 * Pi_AB_inv ;

    K(1:d_heom_A,1:d_heom_A) = -Gamma_L_BA * (G_A_BA_0+G_A_BA_1)  * Gamma_L_AA  - Gamma_R_AB * (G_A_AB_0 + G_A_AB_1) * Gamma_R_AA ;

    % construct BA part
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = ...
        Gamma_R_BA * (G_A_BA_0+G_A_BA_1)  * Gamma_L_AA  + Gamma_L_AB * (G_A_AB_0+G_A_AB_1) * Gamma_R_AA ;

    % K_BB term
    G_B_BA_0 = Pi_BA * (repmat(g_B_BA,[1,1]).* Pi_BA_inv) ;
    G_B_AB_0 = Pi_AB * (repmat(g_B_AB,[1,1]).* Pi_AB_inv) ;
    %     j_B_BA = ((-repmat(g_B_BA,[1,1]) + transpose(repmat(g_B_BA,[1,1])))./(-repmat(Lambda_BA,[1,1])+transpose(repmat(Lambda_BA,[1,1])))) ;
    %     [rows,cols] = find(isnan(j_B_BA)) ;
    %     for r = 1:numel(rows)
    %         j_B_BA(rows(r),cols(r)) = g_B_BA_1(rows(r)) ;
    %     end
    %     %     j_B_BA(isnan(j_B_BA)) = g_B_BA_1(isnan(j_B_BA)) ;
    %     j_B_AB = ((-repmat(g_B_AB,[1,1]) + transpose(repmat(g_B_AB,[1,1])))./(-repmat(Lambda_AB,[1,1])+transpose(repmat(Lambda_AB,[1,1])))) ;
    %     [rows,cols] = find(isnan(j_B_AB)) ;
    %     for r = 1:numel(rows)
    %         j_B_AB(rows(r),cols(r)) = g_B_AB_1(rows(r)) ;
    %     end
    %     j_B_AB(isnan(j_B_AB)) = g_B_AB_1(isnan(j_B_AB));
    %     G_B_BA_1 = Pi_BA*((Pi_BA_inv * L_phonon_BA * Pi_BA).* ((repmat(g_B_BA,[n_ados,1]) - transpose(repmat(g_B_BA,[n_ados,1])))./(-repmat(Lambda_BA,[n_ados,1])+transpose(repmat(Lambda_BA,[n_ados,1])))))*Pi_BA_inv ;
    %     G_B_AB_1 = Pi_AB*((Pi_AB_inv * L_phonon_AB * Pi_AB).* ((repmat(g_B_AB,[n_ados,1]) - transpose(repmat(g_B_AB,[n_ados,1])))./(-repmat(Lambda_AB,[n_ados,1])+transpose(repmat(Lambda_AB,[n_ados,1])))))*Pi_AB_inv ;
    %     G_B_BA_1 = L_phonon_BA *g_B_BA_1 ;
    %     G_B_AB_1 = L_phonon_AB *g_B_AB_1 ;
    %     G_B_BA_1 = Pi_BA*((Pi_BA_inv * L_phonon_BA * Pi_BA).* j_B_BA)*Pi_BA_inv ;
    %     G_B_AB_1 = Pi_AB*((Pi_AB_inv * L_phonon_AB * Pi_AB).* j_B_AB)*Pi_AB_inv ;

    G_B_BA_1 = L_tilde_phonon_BA ;
    [rows,cols] = find(G_B_BA_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            G_B_BA_1(rows(r),cols(r)) = G_B_BA_1(rows(r),cols(r)) * g_B_BA_1(n) ;
        else
            G_B_BA_1(rows(r),cols(r)) = G_B_BA_1(rows(r),cols(r)) * (g_B_BA(n)-g_B_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
        end
    end
    G_B_BA_1 = Pi_BA * G_B_BA_1 * Pi_BA_inv ;

    G_B_AB_1 = L_tilde_phonon_AB ;
    [rows,cols] = find(G_B_AB_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            G_B_AB_1(rows(r),cols(r)) = G_B_AB_1(rows(r),cols(r)) * g_B_AB_1(n) ;
        else
            G_B_AB_1(rows(r),cols(r)) = G_B_AB_1(rows(r),cols(r)) * (g_B_AB(n)-g_B_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
        end
    end
    G_B_AB_1 = Pi_AB * G_B_AB_1 * Pi_AB_inv ;

    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        -Gamma_L_AB * (G_B_AB_0+G_B_AB_1) * Gamma_L_BB  - Gamma_R_BA * (G_B_BA_0+G_B_BA_1) * Gamma_R_BB;

    % construct the AB part
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        Gamma_R_AB * (G_B_AB_0+G_B_AB_1) * Gamma_L_BB  + Gamma_L_BA * (G_B_BA_0+G_B_BA_1) * Gamma_R_BB ;
elseif (AB_coupling_info.method == "second-order phonon NZ" || AB_coupling_info.method == "second-order phonon NZ 2")
    % construct the coherence blocks of the heom liouvillian
    L_phonon_AB = constructHEOMABGenerator(0*H_sys_A,0*H_sys_B,V_As,V_Bs,heom_bath_info,heom_truncation_info,heom_structure_A) ;
    L_phonon_BA = constructHEOMABGenerator(0*H_sys_B,0*H_sys_A,V_Bs,V_As,heom_bath_info,heom_truncation_info,heom_structure_A) ;
    if (AB_coupling_info.method == "second-order phonon NZ 2")
        L_phonon_AB = L_phonon_AB + kron(spdiags([heom_structure_A.ado_gammas],[0],n_ados,n_ados),speye(d_hilb_A*d_hilb_B)) ;
        L_phonon_BA = L_phonon_BA + kron(spdiags([heom_structure_A.ado_gammas],[0],n_ados,n_ados),speye(d_hilb_A*d_hilb_B)) ;

    end

    % construct Gamma operators
    Gamma = AB_coupling_info.coupling_matrix ;
    Gamma_L_AA = kron(id_ados,kron(Gamma',id_hilb_A)) ;
    Gamma_L_BA = kron(id_ados,kron(Gamma,id_hilb_A)) ;
    Gamma_R_AA = kron(id_ados,kron(id_hilb_A,transpose(Gamma))) ;
    Gamma_R_AB = kron(id_ados,kron(id_hilb_A,transpose(Gamma'))) ;
    Gamma_L_BB = kron(id_ados,kron(Gamma,id_hilb_B)) ;
    Gamma_L_AB = kron(id_ados,kron(Gamma',id_hilb_B)) ;
    Gamma_R_BB = kron(id_ados,kron(id_hilb_B,transpose(Gamma'))) ;
    Gamma_R_BA = kron(id_ados,kron(id_hilb_B,transpose(Gamma))) ;

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
    L_tilde_phonon_AB = Pi_AB_inv * L_phonon_AB * Pi_AB ;
    L_tilde_phonon_BA = Pi_BA_inv * L_phonon_BA * Pi_BA ;
    d_AB = d_hilb_A * d_hilb_B ;
    % construct some fourier transforms
    if (AB_coupling_info.method == "second-order phonon NZ")
        g_A_AB = repmat(trapz(ts,conj(c_A_ts).*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_A_BA = repmat(trapz(ts,c_A_ts.*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_B_AB = repmat(trapz(ts,c_B_ts.*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_B_BA = repmat(trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_A_AB_1 = repmat(trapz(ts,ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_A_BA_1 = repmat(trapz(ts,ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_B_AB_1 = repmat(trapz(ts,ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_B_BA_1 = repmat(trapz(ts,ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_A_AB_2 = repmat(0.5*trapz(ts,ts.*ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_A_BA_2 = repmat(0.5*trapz(ts,ts.*ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        g_B_AB_2 = repmat(0.5*trapz(ts,ts.*ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[n_ados,1]) ;
        g_B_BA_2 = repmat(0.5*trapz(ts,ts.*ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[n_ados,1]) ;
        Lambda_AB = repmat(Lambda_AB,[n_ados,1]) ;
        Lambda_BA = repmat(Lambda_BA,[n_ados,1]) ;
    elseif (AB_coupling_info.method == "second-order phonon NZ 2")
        Lambda_AB = repmat(Lambda_AB,[n_ados,1]) - kron(heom_structure_A.ado_gammas,ones([d_hilb_A *d_hilb_B,1]));
        Lambda_BA = repmat(Lambda_BA,[n_ados,1]) - kron(heom_structure_A.ado_gammas,ones([d_hilb_A *d_hilb_B,1]));
        g_A_AB = repmat(trapz(ts,conj(c_A_ts).*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_A_BA = repmat(trapz(ts,c_A_ts.*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_B_AB = repmat(trapz(ts,c_B_ts.*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_B_BA = repmat(trapz(ts,conj(c_B_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_A_AB_1 = repmat(trapz(ts,ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_A_BA_1 = repmat(trapz(ts,ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_B_AB_1 = repmat(trapz(ts,ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_B_BA_1 = repmat(trapz(ts,ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_A_AB_2 = repmat(0.5*trapz(ts,ts.*ts.*conj(c_A_ts).*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_A_BA_2 = repmat(0.5*trapz(ts,ts.*ts.*(c_A_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
        g_B_AB_2 = repmat(0.5*trapz(ts,ts.*ts.*c_B_ts.*exp(Lambda_AB.*ts),2),[1,1]) ;
        g_B_BA_2 = repmat(0.5*trapz(ts,ts.*ts.*conj(c_B_ts).*exp(Lambda_BA.*ts),2),[1,1]) ;
       
    end

    % K_AA term
    G_A_BA_0 = Pi_BA * (repmat(g_A_BA,[1,1]).* Pi_BA_inv) ;
    G_A_AB_0 = Pi_AB * (repmat(g_A_AB,[1,1]).* Pi_AB_inv) ;
    % first and second order phonon terms
    G_A_BA_1 = L_tilde_phonon_BA ;
    G_A_BA_2 = L_tilde_phonon_BA ;
    L_tilde_phonon_BA_deltaLambda = sparse([],[],[],d_AB* n_ados,d_AB* n_ados) ;
    L_tilde_phonon_BA_invDLambda = sparse([],[],[],d_AB* n_ados,d_AB* n_ados) ;
    [rows,cols] = find(G_A_BA_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            G_A_BA_1(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * g_A_BA_1(n) ;
            G_A_BA_2(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * g_A_BA_2(n) ;
            L_tilde_phonon_BA_deltaLambda(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) ;
        else
            G_A_BA_1(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * (g_A_BA(n)-g_A_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
            G_A_BA_2(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * (g_A_BA(n)-g_A_BA(m)-(Lambda_BA(n)-Lambda_BA(m))*g_A_BA_1(m))/((Lambda_BA(n)-Lambda_BA(m))^2) ;
            L_tilde_phonon_BA_invDLambda(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * (1)/(Lambda_BA(n)-Lambda_BA(m)) ;
        end
    end
    X = L_tilde_phonon_BA * L_tilde_phonon_BA_invDLambda ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * g_A_BA_1(n) ;
       
        else
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * (g_A_BA(n)-g_A_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
        end
    end   
    G_A_BA_2 = Pi_BA * (G_A_BA_2 * L_tilde_phonon_BA_deltaLambda + G_A_BA_1*L_tilde_phonon_BA_invDLambda - X) *Pi_BA_inv ;
    G_A_BA_1 = Pi_BA * G_A_BA_1 * Pi_BA_inv ;

    % The AB terms
    G_A_AB_1 = L_tilde_phonon_AB ;
    G_A_AB_2 = L_tilde_phonon_AB ;
    L_tilde_phonon_AB_deltaLambda = sparse([],[],[],d_AB* n_ados,d_AB* n_ados) ;
    L_tilde_phonon_AB_invDLambda = sparse([],[],[],d_AB* n_ados,d_AB* n_ados) ;
    [rows,cols] = find(G_A_AB_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            G_A_AB_1(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * g_A_AB_1(n) ;
            G_A_AB_2(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * g_A_AB_2(n) ;
            L_tilde_phonon_AB_deltaLambda(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) ;
        else
            G_A_AB_1(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * (g_A_AB(n)-g_A_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
            G_A_AB_2(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * (g_A_AB(n)-g_A_AB(m)-(Lambda_AB(n)-Lambda_AB(m))*g_A_AB_1(m))/((Lambda_AB(n)-Lambda_AB(m))^2) ;
            L_tilde_phonon_BA_invDLambda(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * (1)/(Lambda_AB(n)-Lambda_AB(m)) ;
        end
    end
    X = L_tilde_phonon_AB * L_tilde_phonon_AB_invDLambda ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * g_A_AB_1(n) ;
       
        else
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * (g_A_AB(n)-g_A_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
        end
    end   
    G_A_AB_2 = Pi_AB * (G_A_AB_2 * L_tilde_phonon_AB_deltaLambda + G_A_AB_1*L_tilde_phonon_AB_invDLambda - X) *Pi_AB_inv ;
    G_A_AB_1 = Pi_AB * G_A_AB_1 * Pi_AB_inv ;
%     G_A_AB_1 = L_tilde_phonon_AB ;
%     [rows,cols] = find(G_A_AB_1) ;
%     for r = 1:numel(rows)
%         n = rows(r) ;
%         m = cols(r) ;
%         if (Lambda_AB(n)==Lambda_AB(m))
%             G_A_AB_1(rows(r),cols(r)) = G_A_AB_1(rows(r),cols(r)) * g_A_AB_1(n) ;
%         else
%             G_A_AB_1(rows(r),cols(r)) = G_A_AB_1(rows(r),cols(r)) * (g_A_AB(n)-g_A_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
%         end
%     end
%     G_A_AB_1 = Pi_AB * G_A_AB_1 * Pi_AB_inv ;

    K(1:d_heom_A,1:d_heom_A) = -Gamma_L_BA * (G_A_BA_0+G_A_BA_1 +G_A_BA_2)  * Gamma_L_AA  - Gamma_R_AB * (G_A_AB_0 + G_A_AB_1 + G_A_AB_2) * Gamma_R_AA ;

    % construct BA part
    K((d_heom_A+1):(d_heom_A+d_heom_B),1:d_heom_A) = ...
        Gamma_R_BA * (G_A_BA_0+G_A_BA_1+G_A_BA_2)  * Gamma_L_AA  + Gamma_L_AB * (G_A_AB_0+G_A_AB_1+G_A_AB_2) * Gamma_R_AA ;

    % K_BB term
    G_B_BA_0 = Pi_BA * (repmat(g_B_BA,[1,1]).* Pi_BA_inv) ;
    G_B_AB_0 = Pi_AB * (repmat(g_B_AB,[1,1]).* Pi_AB_inv) ;
    % first and second order phonon terms
    G_B_BA_1 = L_tilde_phonon_BA ;
    G_B_BA_2 = L_tilde_phonon_BA ;
    [rows,cols] = find(G_B_BA_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            G_B_BA_1(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * g_B_BA_1(n) ;
            G_B_BA_2(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * g_B_BA_2(n) ;
        else
            G_B_BA_1(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * (g_B_BA(n)-g_B_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
            G_B_BA_2(rows(r),cols(r)) = L_tilde_phonon_BA(rows(r),cols(r)) * (g_B_BA(n)-g_B_BA(m)-(Lambda_BA(n)-Lambda_BA(m))*g_B_BA_1(m))/((Lambda_BA(n)-Lambda_BA(m))^2) ;
        end
    end
    X = L_tilde_phonon_BA * L_tilde_phonon_BA_invDLambda ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_BA(n)==Lambda_BA(m))
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * g_B_BA_1(n) ;
       
        else
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * (g_B_BA(n)-g_B_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
        end
    end   
    G_B_BA_2 = Pi_BA * (G_B_BA_2 * L_tilde_phonon_BA_deltaLambda + G_B_BA_1*L_tilde_phonon_BA_invDLambda - X) *Pi_BA_inv ;
    G_B_BA_1 = Pi_BA * G_B_BA_1 * Pi_BA_inv ;

    % The AB terms
    G_B_AB_1 = L_tilde_phonon_AB ;
    G_B_AB_2 = L_tilde_phonon_AB ;
    [rows,cols] = find(G_B_AB_1) ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            G_B_AB_1(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * g_B_AB_1(n) ;
            G_B_AB_2(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * g_B_AB_2(n) ;
        else
            G_B_AB_1(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * (g_B_AB(n)-g_B_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
            G_B_AB_2(rows(r),cols(r)) = L_tilde_phonon_AB(rows(r),cols(r)) * (g_B_AB(n)-g_B_AB(m)-(Lambda_AB(n)-Lambda_AB(m))*g_B_AB_1(m))/((Lambda_AB(n)-Lambda_AB(m))^2) ;
        end
    end
    X = L_tilde_phonon_AB * L_tilde_phonon_AB_invDLambda ;
    for r = 1:numel(rows)
        n = rows(r) ;
        m = cols(r) ;
        if (Lambda_AB(n)==Lambda_AB(m))
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * g_B_AB_1(n) ;
       
        else
            X(rows(r),cols(r)) = X(rows(r),cols(r)) * (g_B_AB(n)-g_B_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
        end
    end   
    G_B_AB_2 = Pi_AB * (G_B_AB_2 * L_tilde_phonon_AB_deltaLambda + G_B_AB_1*L_tilde_phonon_AB_invDLambda - X) *Pi_AB_inv ;
    G_B_AB_1 = Pi_AB * G_B_AB_1 * Pi_AB_inv ;
%       G_B_BA_1 = L_tilde_phonon_BA ;
%     [rows,cols] = find(G_B_BA_1) ;
%     for r = 1:numel(rows)
%         n = rows(r) ;
%         m = cols(r) ;
%         if (Lambda_BA(n)==Lambda_BA(m))
%             G_B_BA_1(rows(r),cols(r)) = G_B_BA_1(rows(r),cols(r)) * g_B_BA_1(n) ;
%         else
%             G_B_BA_1(rows(r),cols(r)) = G_B_BA_1(rows(r),cols(r)) * (g_B_BA(n)-g_B_BA(m))/(Lambda_BA(n)-Lambda_BA(m)) ;
%         end
%     end
%     G_B_BA_1 = Pi_BA * G_B_BA_1 * Pi_BA_inv ;
% 
%     G_B_AB_1 = L_tilde_phonon_AB ;
%     [rows,cols] = find(G_B_AB_1) ;
%     for r = 1:numel(rows)
%         n = rows(r) ;
%         m = cols(r) ;
%         if (Lambda_AB(n)==Lambda_AB(m))
%             G_B_AB_1(rows(r),cols(r)) = G_B_AB_1(rows(r),cols(r)) * g_B_AB_1(n) ;
%         else
%             G_B_AB_1(rows(r),cols(r)) = G_B_AB_1(rows(r),cols(r)) * (g_B_AB(n)-g_B_AB(m))/(Lambda_AB(n)-Lambda_AB(m)) ;
%         end
%     end
%     G_B_AB_1 = Pi_AB * G_B_AB_1 * Pi_AB_inv ;

    K((d_heom_A+1):(d_heom_A+d_heom_B),(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        -Gamma_L_AB * (G_B_AB_0+G_B_AB_1+G_B_AB_2) * Gamma_L_BB  - Gamma_R_BA * (G_B_BA_0+G_B_BA_1+G_B_BA_2) * Gamma_R_BB;

    % construct the AB part
    K(1:d_heom_A,(d_heom_A+1):(d_heom_A+d_heom_B)) = ...
        Gamma_R_AB * (G_B_AB_0+G_B_AB_1+G_B_AB_2) * Gamma_L_BB  + Gamma_L_BA * (G_B_BA_0+G_B_BA_1+G_B_BA_2) * Gamma_R_BB ;


end

if AB_coupling_info.phonon_method == "simplified"
    % first the
end


end