function [H,H_no_ls,H_R,H_q,R,q] = constructCavityHamiltonian(omega_cav,M,chi,V_params,mu_params,n_R,n_q,R_range,q_range,eta,omega_D)

% grids for q and R
R = ((0:(n_R-1))/(n_R-1) * (R_range(2)-R_range(1)) + R_range(1))' ;
q = ((0:(n_q-1))/(n_q-1) * (q_range(2)-q_range(1)) + q_range(1))' ;
dR = R(2)-R(1) ;
dq = q(2)-q(1) ;
d = n_R * n_q ;


% V_0(R)
V_0_R = V_params(1)*(R.^2) + V_params(2) * (R.^4) + V_params(3)  + V_params(4) * R ;

% mu(R)
mu_R = mu_params(1) * tanh(mu_params(2)*R) + mu_params(3) * R ;

% V(R,q)
V = kron(V_0_R,ones([n_q,1]))  ...
    + (omega_cav^2/2) * kron(ones([n_R,1]),q.^2) ...
    + (omega_cav^2/2) * (2/(omega_cav^3))*(chi^2) * kron(mu_R.*mu_R,ones([n_q,1])) ...
    + (omega_cav^2) * sqrt((2/(omega_cav^3))) * chi * kron(mu_R,q);
V_op = spdiags([V],0,d,d) ;

% T 2nd order
T_R = (-0.5/(M*dR*dR)) * spdiags([1,-2,1].*ones([n_R,1]),-1:1,n_R,n_R) ;
T_q = (-omega_cav^0 / (2*dq*dq)) * spdiags([1,-2,1].*ones([n_q,1]),-1:1,n_q,n_q) ;
% 4th order finite diff
T_R = (-0.5/(M*dR*dR)) * spdiags([-1/12,4/3,-5/2,4/3,-1/12].*ones([n_R,1]),-2:2,n_R,n_R) ;
T_q = (-omega_cav^0 / (2*dq*dq)) * spdiags([-1/12,4/3,-5/2,4/3,-1/12].*ones([n_q,1]),-2:2,n_q,n_q) ;
% Colbert-Miller DVR
T_R = (-0.5/(M)) *  constructCMDVR1DSecondDerivative(dR,n_R) ; 
T_q = (-omega_cav^0 / (2)) * constructCMDVR1DSecondDerivative(dq,n_q) ; 
T_op = kron(T_R,speye(n_q)) + kron(speye(n_R),T_q) ;

H = T_op + V_op ;

V_ls = (eta * omega_D/(M*pi)) * kron(R.^2,ones([n_q,1])) ;
H_ls_op = spdiags([V_ls],0,d,d) ; 
H_no_ls = H ;
H = H_no_ls + H_ls_op ;

H_R = spdiags(V_0_R,0,n_R,n_R) + T_R ;
H_q = spdiags((omega_cav^2/2) *q.^2,0,n_q,n_q) + T_q ;

end