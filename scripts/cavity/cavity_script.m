% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters
omega_cav = 0.006269431 ;
M = 1836.0 ;
chi = 1*0.00234562 ;
mu_params = [-1.90249,1.26426,0.37044] ;
V_params = [-0.021087856, 0.0033107783,0.033160555,3.6749309e-6] ;
% bath parameters
beta = 1052.584412992859 ;
beta = 100 ;
% debye bath parameters
eta =  6.601876175e-8 ; 
omega_D = 0.006269431 ;
% omega_D = 0.01 ;
lambda_D = eta / (2 * omega_D) ;

n_q = 70 ;
n_R = 100 ;
d = n_q * n_R ;
R_range = [-5,5] ;
q_range = [-100,100] ;
[H,H_no_ls,H_R,H_q,R,q] = constructCavityHamiltonian(omega_cav,M,chi,V_params,mu_params,n_R,n_q,R_range,q_range,eta,omega_D) ;

n_E = 50 ;
[Psi_E,E_mat] = eigs(H , n_E, 'smallestabs') ;
E = diag(E_mat) - E_mat(1,1) ;

H_sys_E = sparse(E_mat) ;
V_E = Psi_E' * spdiags(kron(R,ones([n_q,1])),0,n_q*n_R,n_q*n_R) * Psi_E ; 
H_sys  = H_sys_E ;
V = V_E ;
% [Psi_DVR,V_DVR] = eig(V_E) ;
% H_sys = Psi_DVR' * H_sys_E * Psi_DVR ;
% V = sparse(V_DVR) ;

% rho_0
n_init = 10 ;
rho_0_sys = zeros([n_E,n_E]) ;
rho_0_sys(n_init,n_init) = 1.0  ;
% rho_0_sys = Psi_DVR' * rho_0_sys * Psi_DVR ;

% O_sys 
P_1 = sparse(zeros([n_E,n_E]));
P_1(1,1) = 1.0 ;
P_2 = sparse(zeros([n_E,n_E]));
P_2(2,2) = 1.0 ;
% P_1 = Psi_DVR' * P_1 * Psi_DVR ;
% P_2 = Psi_DVR' * P_2 * Psi_DVR ;
O_sys = {H_sys,P_1,P_2} ;

% HEOM truncation
L_max = 1 ;
M_max = 0 ;

% dynamics
dt = 10.0 ;
n_steps = 100000 ;
krylov_dim = 16 ;
krylov_tol = 1e-10 ;
low_temp_corr_method = "low temperature correction" ; % NZ2 method
low_temp_corr_method = "NZ2" ;

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = H_sys ;
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",V,...
    "spectral_density","debye (pade)","omega_D",omega_D,"lambda_D",lambda_D,...
    "N_pade",M_max,"approximant_type","[N/N]")} ;
% full_system.baths = {struct("V",V,...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.beta = beta ;

% a struct that contains information about the HEOM dynamics
heom_dynamics = struct() ;
% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct ;
heom_dynamics.integrator.method = "SIA" ;
heom_dynamics.integrator.dt = dt ;
heom_dynamics.integrator.n_steps = n_steps ;
heom_dynamics.integrator.krylov_dim = krylov_dim ;
heom_dynamics.integrator.krylov_tol = krylov_tol ;

% hierarchy trunction information
heom_dynamics.heom_truncation = struct() ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
heom_dynamics.heom_truncation.heom_termination = low_temp_corr_method ;
heom_dynamics.heom_truncation.diagonal_only_term = true ;
heom_dynamics.heom_truncation.termination_k_max = 200 ;


% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
n_t_plot = min([500,n_steps]) ;
skip = floor(n_steps/n_t_plot);
O_t = O_t(:,1:skip:end) ;
t = t(1:skip:end) ;


% plot population as a function of time
% figure 
% plot(t,0.5+0.5*O_t(3,:))