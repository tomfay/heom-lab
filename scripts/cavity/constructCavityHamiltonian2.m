function [H,H_no_ls,H_R,H_q,R] = constructCavityHamiltonian2(omega_cav,M,chi,V_params,mu_params,n_R,n_q,R_range,eta,omega_D)

% grids for q and R
R = ((0:(n_R-1))/(n_R-1) * (R_range(2)-R_range(1)) + R_range(1))' ;
% q = ((0:(n_q-1))/(n_q-1) * (q_range(2)-q_range(1)) + q_range(1))' ;
dR = R(2)-R(1) ;
% dq = q(2)-q(1) ;
d = n_R * n_q ;


% V_0(R)
V_0_R = V_params(1)*(R.^2) + V_params(2) * (R.^4) + V_params(3)  + V_params(4) * R ;

% mu(R)
mu_R = mu_params(1) * tanh(mu_params(2)*R) + mu_params(3) * R ;
mu_R_op = spdiags([mu_R],[0],n_R,n_R) ;

% H_cav 
H_q_0 = spdiags([omega_cav*(0.5+(0:(n_q-1))')],[0],n_q,n_q) ;

% q
q_op = spdiags((sqrt(0.5/omega_cav))*[sqrt(1:n_q)',sqrt(0:(n_q-1))'],[-1,1],n_q,n_q) ;

% V(R,q)
V = kron(V_0_R,ones([n_q,1]))  ...
    + (omega_cav^2/2) * (2/(omega_cav^3))*(chi^2) * kron(mu_R.*mu_R,ones([n_q,1])) ;
    
V_op = spdiags([V],0,d,d) + (omega_cav^2) * sqrt((2/(omega_cav^3))) * chi * kron(mu_R_op,q_op) +  kron(speye(n_R),H_q_0);

% T 2nd order
T_R = (-0.5/(M*dR*dR)) * spdiags([1,-2,1].*ones([n_R,1]),-1:1,n_R,n_R) ;
% T_q = (-omega_cav^0 / (2*dq*dq)) * spdiags([1,-2,1].*ones([n_q,1]),-1:1,n_q,n_q) ;
% 4th order finite diff
T_R = (-0.5/(M*dR*dR)) * spdiags([-1/12,4/3,-5/2,4/3,-1/12].*ones([n_R,1]),-2:2,n_R,n_R) ;
% T_q = (-omega_cav^0 / (2*dq*dq)) * spdiags([-1/12,4/3,-5/2,4/3,-1/12].*ones([n_q,1]),-2:2,n_q,n_q) ;
% Colbert-Miller DVR
T_R = (-0.5/(M)) *  constructCMDVR1DSecondDerivative(dR,n_R) ; 
% T_q = (-omega_cav^0 / (2)) * constructCMDVR1DSecondDerivative(dq,n_q) ; 
T_op = kron(T_R,speye(n_q)) ;

H = T_op + V_op ;

V_ls = (eta * omega_D/(M*pi)) * kron(R.^2,ones([n_q,1])) ;
H_ls_op = spdiags([V_ls],0,d,d) ; 
H_no_ls = H ;
H = H_no_ls + H_ls_op ;

H_R = spdiags(V_0_R,0,n_R,n_R) + T_R ;
H_q = H_q_0 ;

end