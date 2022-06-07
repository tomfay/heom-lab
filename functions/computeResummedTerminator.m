function T = computeResummedTerminator(V_j_L,V_j_R,c_jk,cbar_jk,nu_jk,n_jk,gamma_ado,n_max,Pi_sys,Pi_sys_inv,lambda_sys,L_sys)
% computes a resummed contribution to the terminator

L_minus =  ((-1.0i*c_jk/sqrt(abs(c_jk))) * V_j_L - (-1.0i*conj(cbar_jk)/sqrt(abs(c_jk))) * V_j_R) ;
L_plus = (-1.0i*sqrt(abs(c_jk))) * (V_j_L - V_j_R) ;

L_plus_inv = (pinv(full(L_plus))) ;
L_minus_inv = (pinv(full(L_minus))) ;

d = size(Pi_sys,1) ;

% The K term
K_m_minus_1 = zeros([d,d]) ; % K_0
K_m = eye(d) ;                  % K_1 
K_m_plus_1 = zeros([d,d]) ;

% compute K_{n_max+1}
m = 1 ;
C_m_inv = (1/sqrt(n_jk + m)) * L_plus_inv * (L_sys - speye(d)*(gamma_ado + nu_jk * m)) ;
A_m = eye(d) ;
K_m_plus_1 = - C_m_inv * (-K_m + A_m * K_m_minus_1) ; % K_2
for m = 2:(n_max-1) 
    K_m_minus_1 = K_m ;
    K_m = K_m_plus_1 ;
    
    C_m_inv = (1/sqrt(n_jk +  m)) * L_plus_inv * (L_sys - speye(d)*(gamma_ado + nu_jk * m)) ;
    G_m = Pi_sys * ( (1./(-lambda_sys + (gamma_ado + nu_jk * m))).*Pi_sys_inv) ;
    A_m = sqrt(n_jk+m) * (G_m * L_minus) ;
    K_m_plus_1 = - C_m_inv * (-K_m + A_m * K_m_minus_1) ;
    
end
m = n_max ;
K_m_minus_1 = K_m ;
K_m = K_m_plus_1 ;
C_m_inv = eye(d) ;
G_m = Pi_sys * ( (1./(-lambda_sys + (gamma_ado + nu_jk * m))).*Pi_sys_inv) ;
A_m = sqrt(n_jk+m) * (G_m * L_minus) ;
K_m_plus_1 = - C_m_inv * (-K_m + A_m * K_m_minus_1) ;


% compute N_1 
N_m_minus_1 = zeros([d,d]) ;
N_m = eye(d) ;
N_m_plus_1 = zeros([d,d]) ;
m = n_max ;
A_m_plus_1 = eye(d) ;
C_m_minus_1_inv = (1/sqrt(n_jk + m )) * L_plus_inv * (L_sys - speye(d)*(gamma_ado + nu_jk * (m-1))) ;
N_m_minus_1 = - (N_m  + N_m_plus_1 * A_m_plus_1) *C_m_minus_1_inv ;
for m = (n_max-1):-1:2 
    N_m_plus_1 = N_m ;
    N_m = N_m_minus_1 ;
    G_m_plus_1 = Pi_sys * ( (1./(-lambda_sys + (gamma_ado + nu_jk * (m+1)))).*Pi_sys_inv) ;
    A_m_plus_1 = sqrt(n_jk+m+1) * (G_m_plus_1 * L_minus) ;
    C_m_minus_1_inv = (1/sqrt(n_jk + m  )) * L_plus_inv * (L_sys - speye(d)*(gamma_ado + nu_jk * (m-1))) ;
    N_m_minus_1 = - (N_m  + N_m_plus_1 * A_m_plus_1) *C_m_minus_1_inv ;
end

G_eff = pinv(K_m_plus_1) * N_m_minus_1 ;
if(n_max == 1)
    G_eff = eye(d) ;
end
% G_eff
G_m = Pi_sys * ( (1./(-lambda_sys + (gamma_ado + nu_jk))).*Pi_sys_inv)  ;
% G_m = -inv(L_sys - (gamma_ado+nu_jk)*eye(d))
T = (n_jk+1)*L_plus * G_eff * G_m * L_minus ;

end