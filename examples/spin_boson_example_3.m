% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters
epsilon = 0.5 ;
Delta = 1.0 ;
% bath parameters
beta = 0.25 ;
% debye bath parameters
lambda_D = 0.06*Delta/(2*pi) ;
omega_D = 0.05*Delta ;
alpha = pi/4 ;
phi = -pi/2 ;
theta = 3*pi/8 + pi;

% dynamics information
dt = 5e-2 ;
n_steps = 2400 ;
krylov_dim = 8 ;
krylov_tol = 1e-8 ;
L_max = 10 ;
M_max = 1 ;
% Gamma_cut = 4.1*omega_D ;

% matrices of system observable operators to be returned, sigma_x, sigma_y
% sigma_z, and 1
O_sys = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;

% initial state of the system
psi = [cos(alpha) ; sin(alpha)*exp(-1.0i*phi)] ;
rho_0_sys = psi*(psi') ;

% two objects are supplied to the HEOM dynamics function:
% "full_system" specifies the full Hamiltonian (system + bath) and the
% temperature.
% "heom_dynamics" specifies the HEOM truncation, integrator for the
% dynamics and the total propagation time and observables to be calculated.

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon/2,Delta/2];
                     [Delta/2,-epsilon/2]];
% baths is a cell array of structs describign each bath
S_op = [[sin(theta),cos(theta)];[cos(theta),-sin(theta)]] ;
full_system.baths = {struct("V",S_op,...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.beta = beta ;

% a struct that contains information about the HEOM dynamics
heom_dynamics = struct ;
% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct ;
heom_dynamics.integrator.method = "adaptive SIA" ;
heom_dynamics.integrator.dt = dt ;
heom_dynamics.integrator.n_steps = n_steps ;
heom_dynamics.integrator.krylov_dim = krylov_dim ;
heom_dynamics.integrator.krylov_tol = krylov_tol ;

% hierarchy trunction information
heom_dynamics.heom_truncation = struct ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
heom_dynamics.heom_truncation.heom_termination = "NZ2" ;
heom_dynamics.heom_truncation.termination_k_max = 100 ;
heom_dynamics.heom_truncation.diagonal_only_term = true ;

% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
c_0 = 2.99792458e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;