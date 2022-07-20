% A test script for the spin boson model
% In this example the dynamics for a spin boson model with a Debye bath are
% calculated

% Parameters for the problem
% system hamiltonian parameters
epsilon = 14780 ;
epsilon = 14980 ;
% epsilon = 10 ;
% Delta = 5 ;
epsilon = 30 ;
% bath parameters
% beta = 1.0/208.50907518 ;
beta = 1.0 ;
% debye bath parameters
% lambda_D = 220 ;
lambda_D = 0.15 ;
% omega_D = 353.6777 ;
omega_D = 1.0 ;
kappa = 0.2 ;
% lambda_D = 37 ;
% omega_D = 30 ;
Omega_B = 7.5 ;
gamma_B = 1 ;
% lambda_B = 0.5 ;

% dynamics information
dt = 1e-1 ;
n_steps = 200000 ;
krylov_dim = 16 ;
krylov_tol = 1e-10 ;
L_max = 2 ;
M_max = 2 ;
low_temp_corr_method = "low temperature correction" ; % Ishizaki-Tanimura
% low_temp_corr_method = "ow temperature correction NZ2" ; % NZ2 method
% Gamma_cut = (2*pi/beta)*4.01 ;
% Gamma_cut = Omega_B * 2.01 ;
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
full_system.H_sys = [[0,0];
                     [0,-epsilon]];
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",[[0,0];[0,1]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D),...
%     struct("V",[[0,1];[1,0]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",kappa*lambda_D)} ;
struct("V",[[0,1];[1,0]],...
    "spectral_density","UBO","Omega",Omega_B,"lambda",kappa*lambda_D,...
    "gamma",gamma_B)} ;
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
heom_dynamics.heom_truncation = struct() ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
heom_dynamics.heom_truncation.heom_termination = low_temp_corr_method ;
% heom_dynamics.heom_truncation.heom_termination = "low temp correction NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "low temp correction" ;
% heom_dynamics.heom_truncation.heom_termination = "NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "RF2" ;
% heom_dynamics.heom_truncation.heom_termination = "partial resummed" ;
% heom_dynamics.heom_truncation.n_max_resum = 1 ;
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
c_0 = 2.99792458e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;

% plot population as a function of time
figure 
plot(t,0.5+0.5*O_t(3,:))