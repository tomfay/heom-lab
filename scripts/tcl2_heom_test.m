% This script tests the TCL2 approach for integrating out a strongly
% coupled bath from an exciton dimer coupled to a CT state

% Here the |A,k> states are |D*A> and |DA*> and the |B> state is |CT>.
% The |A,k> states are each coupled weakly to a debye bath.
% The |B> state is coupled strongly to a debye bath.

% set up system parameters
% H_s,A parameters
delta_epsilon = 1.0e0 ;
J = 1.0e0 ; 

% set up explicit bath parameters
lambda_D = 1e-1 ;
omega_D = 0.25 ;
beta = 1.0 ;
lambda_AB = 6.0 ;
omega_AB = 0.25 ;
Gamma_AB = 0.1 ;
Delta_E_AB = 5 ;

% dynamics information
dt = 0.01e0 ;
n_steps = 400000 ;
krylov_dim = 8 ;
krylov_tol = 1e-11 ;
Gamma_cut = 16.0 ;
p = 2 ;
L_cut = 16 ;

% parameters for evaluating to AB correlation function
t_max = sqrt((beta/lambda_AB)*log(1/1e-7)) ;
n_t = 1000 ;
n_modes = 512 ; % number of modes used to discretise the spectral density

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys_A = [[delta_epsilon/2,J];
                     [J,-delta_epsilon/2]];
full_system.H_sys_B = [[0]] ;
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V_A",[[1,0];[0,0]],"V_B",[[0]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.baths = [full_system.baths,...
    {struct("V_A",[[0,0];[0,1]],"V_B",[[0]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)}] ;
% full_system.baths = [full_system.baths,...
%     {struct("V_A",[[1,0.5];[0.5,0]],"V_B",[[0]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",0.25*lambda_D)}] ;
% full_system.baths = {struct("V_A",[[1,0];[0,-1]],"V_B",[[0]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.beta = beta ;

% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition
heom_dynamics.rho_0_sys_A = [[1,0];[0,0]] ;
heom_dynamics.rho_0_sys_B = [[0]] ;

% set up observable arrays
heom_dynamics.observables = struct() ;
heom_dynamics.observables.system_A = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;
heom_dynamics.observables.system_B = {[[1]]} ;

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
heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
heom_dynamics.heom_truncation.heom_termination = "markovian" ;
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "coupling weighted cut-off" ;
% heom_dynamics.heom_truncation.L_cut = L_cut ;
% heom_dynamics.heom_truncation.p = p ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;

% details of the strongly coupled bath
AB_coupling_info = struct() ;
AB_coupling_info.baths = {struct("spectral_density","debye","omega_D",omega_AB,"lambda_D",lambda_AB)} ;
AB_coupling_info.beta = beta ;
AB_coupling_info.coupling_matrix = [0;Gamma_AB] ; % Gamma_a,b matrix
AB_coupling_info.Delta_E_AB = Delta_E_AB ;
AB_coupling_info.t_max = t_max ;
AB_coupling_info.n_t = n_t ;
AB_coupling_info.n_modes = n_modes ;
% AB_coupling_info.method = "simplified" ;
% AB_coupling_info.method = "include H_sys" ;
AB_coupling_info.method = "include H_sys NZ" ;
% AB_coupling_info.method = "full NZ" ;
% AB_coupling_info.method = "first-order phonon NZ" ;
% AB_coupling_info.method = "first-order phonon NZ 2" ;
AB_coupling_info.phonon_method = "none" ;

% run the dynamics
[O_t_AB,t_AB,junk] = runHEOMTC2ABDynamics(full_system,heom_dynamics,AB_coupling_info) ;