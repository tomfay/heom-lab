% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters

% potential parameters
% the potential has the form:
% V(q) = (V_0/q_0^4)(q_q_0)^2(q+q_0)^2 + (Delta_E/(2 q_0)) q
planck_const = 6.62607015e-34 ; % J s
planck_const_eV = 4.135667696e-15 ;
speed_of_light = 299792458e2 ; % cm s^-1
mu = 1.67262192369e-27 ;
J_per_eV = 1.602176634e-19 ;
V_0_eV = 1.0 ;
V_0 = V_0_eV / (planck_const_eV * speed_of_light)  ; % barrier in absence of bias in cm^-1
% V_0 = 10000 ;
omega_0 = 3500 ; % well frequency in abscence of bias in cm^-1
Delta_E_eV = 0.25 ; % approximate difference betwen well minima

omega_0_angfreq = 2*pi*speed_of_light * omega_0 ;
k_0_eV = mu * omega_0_angfreq^2  / J_per_eV ;
q_0 = 1e10*sqrt(8*V_0_eV/k_0_eV) ; % q_0 in Angstrom

% DVR grid parameters 
q_min = -2*q_0 ; q_max = 2.0*q_0 ;
n_q = 500 ;
q = ((0:(n_q-1))/(n_q-1) * (q_max-q_min) + q_min)' ;
dq = q(2) - q(1) ;

% construct the Hamiltonian
% kinetic energy
T_op = (1/J_per_eV)*(-planck_const^2/((2*pi)^2*2*mu)) * constructCMDVR1DSecondDerivative(dq*1e-10,n_q) ;
% potential energy
V_q = V_0_eV / (q_0^4) * (q-q_0).^2 .* (q+q_0).^2 + (Delta_E_eV/(2*q_0))*q;
V_op = diag(V_q) ;
H_op = T_op + V_op ;
H_op = (1/(planck_const_eV*speed_of_light)) * H_op ; % convert to cm^-1
% n_E = 10 ;
% [Psi_E,E_mat] = eigs(H_op , n_E, 'smallestabs') ;
% E = diag(E_mat) ;


% debye bath parameters 
beta = 1.0/200.0 ;
omega_D = 100.0 ;
lambda_D = 10.0 ;

% truncation of the system hilbert space to the lowest energy eigenstates
n_E = 10 ;
[Psi_E,E_mat] = eigs(H_op , n_E, 'smallestabs') ;
E = diag(E_mat) - E_mat(1,1) ;

% conversion of the system Hamiltonian to a DVR basis by diagonalising the
% system-bath coupling operator
H_sys_E = sparse(E_mat) ;
V_sys_bath = spdiags((q/q_0),0,n_q,n_q) ;
V_E = Psi_E' * V_sys_bath * Psi_E ; 
H_sys  = H_sys_E ;
V = V_E ;
[Psi_DVR,V_DVR] = eig(V_E) ;
H_sys = Psi_DVR' * H_sys_E * Psi_DVR ;
V = sparse(V_DVR) ;

% rho_0 - system is initialised in the second energy eigenstate, localised
% in the higher energy well
n_init = 2 ;
rho_0_sys = zeros([n_E,n_E]) ;
rho_0_sys(n_init,n_init) = 1.0  ;
rho_0_sys = Psi_DVR' * rho_0_sys * Psi_DVR ;

% O_sys system observables
O_sys = {} ;
for n = 1:n_E 
P_n = sparse(zeros([n_E,n_E]));
P_n(n,n) = 1 ;
P_n = Psi_DVR' * P_n * Psi_DVR ;
O_sys = [O_sys,{P_n}] ;
end
O_sys = [O_sys,{H_sys}] ;


% HEOM truncation
L_max = 3 ;
M_max = 2 ;

% dynamics parameters 
dt = 1e-5 ;
n_steps = 4000 ;
krylov_dim = 16 ;
krylov_tol = 1e-11 ;
low_temp_corr_method = "low temp correction NZ2" ; % NZ2 method
% low_temp_corr_method = "NZ2" ;

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = H_sys ;
% baths is a cell array of structs describign each bath
% full_system.baths = {struct("V",V,...
%     "spectral_density","debye (pade)","omega_D",omega_D,"lambda_D",lambda_D,...
%     "N_pade",M_max,"approximant_type","[N/N]")} ;
full_system.baths = {struct("V",V,...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
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
heom_dynamics.observables = struct() ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
n_t_plot = min([500,n_steps]) ;
skip = floor(n_steps/n_t_plot);
O_t = O_t(:,1:skip:end) ;
t = t(1:skip:end) ;
% time in femtoseconds
t_fs = 1e12*t / (2*pi*speed_of_light) ;


% plot population as a function of time
% figure 
% plot(t,0.5+0.5*O_t(3,:))