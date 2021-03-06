% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters
% change as appropriate
epsilon = 20 ;
Delta = 5 ;
% bath parameters
beta = 1.0 ;
% debye bath parameters
lambda_D = 0.1 ;
omega_D = 1.0 ;

% set the choice of terminator Choose as appropriate
terminator = "low temp correction" ; % original Ishizaki-Tanimura correction
terminator = "NZ2" ; % Fay's Zwanzig projetion based correction
% terminator = "low temp correction NZ2" ;
% terminator = "none" ;

% dynamics information
dt = 1e-2 ;
n_steps = 200000 ;
krylov_dim = 16 ;
krylov_tol = 1e-10 ;
L_max = 3 ;
M_max = 2 ;
Gamma_cut = 6.01*pi/beta ;
Gamma_cut = 80*omega_D ;
% matrices of system observable operators to be returned, sigma_x, sigma_y
% sigma_z, and 1
O_sys = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;

% initial state of the system
rho_0_sys = [[1,0];[0,0]] ;

% two objects are supplied to the HEOM dynamics function:
% "full_system" specifies the full Hamiltonian (system + bath) and the
% temperature.
% "heom_dynamics" specifies the HEOM truncation, integrator for the
% dynamics and the total propagation time and observables to be calculated.

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon/2,Delta];
                     [Delta,-epsilon/2]];
% baths is a cell array of structs describign each bath
% full_system.baths = {struct("V",[[1,0];[0,-1]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.baths = {struct("V",[[1,0];[0,-1]],...
    "spectral_density","debye (pade)","omega_D",omega_D,"lambda_D",lambda_D,"N_pade",M_max,"approximant_type","[N/N]")} ;
% full_system.baths = {struct("V",[[1,0];[0,0]],...
%     "spectral_density","UBO","Omega",Omega_B,"lambda",lambda_B,...
%     "gamma",gamma_B)} ;
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
heom_dynamics.heom_truncation = struct ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
% heom_dynamics.heom_truncation.heom_termination = "low temp correction NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "low temp correction" ;
% heom_dynamics.heom_truncation.heom_termination = "NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "RF2" ;
% heom_dynamics.heom_truncation.n_max_resum = 1 ;
heom_dynamics.heom_truncation.heom_termination = terminator ;
heom_dynamics.heom_truncation.diagonal_only_term = true ;
heom_dynamics.heom_truncation.termination_k_max = 200 ;


% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
n_t_plot = min([5000,n_steps]) ;
skip = floor(n_steps/n_t_plot);
O_t_full = O_t ;
t_full = t ;
O_t = O_t(:,1:skip:end) ;
t = t(1:skip:end) ;
c_0 = 2.99792458e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;