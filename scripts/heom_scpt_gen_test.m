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
lambda_D = 1.0e0 ;
omega_D = 1.75 ;
beta = 1.0 ;
lambda_12 = 5.0 ;
omega_12 = 1.75 ;
lambda_13 = 5.0 ;
omega_13 = 1.75 ;
Gamma_12 = 0.5 ;
Gamma_13 = 0.0 ;
Gamma_23 = 0.5 ;
omega_23 = 1.75 ;
lambda_23 = lambda_12 ;
E_1 = 0 ; E_2 = -2 ; E_3 = -10 ;

% dynamics information
dt = 0.01e0 ;
n_steps = 100000 ;
krylov_dim = 9 ;
krylov_tol = 1e-12 ;
Gamma_cut = 10*omega_D ;
Gamma_cut_trunc = 0.75 * Gamma_cut ;
p = 1 ;
L_cut = 5.0 ;

% parameters for evaluating to AB correlation fction
t_max = sqrt((beta/lambda_12)*log(1/1e-9)) ;
n_t = 400 ;
n_modes = 512 ; % number of modes used to discretise the spectral density

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = {} ;
full_system.H_sys{1} = [[delta_epsilon/2,J];
                     [J,-delta_epsilon/2]];
full_system.H_sys{2} = [[0]] ;
full_system.H_sys{3} = [[0]] ;

% baths is a cell array of structs describign each bath
full_system.baths = {} ;
full_system.baths = {struct(...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.baths = [full_system.baths,...
    {struct(...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)}] ;
full_system.Vs = {} ;
full_system.Vs{1} = {[[1,0];[0,0]],[[1]],[[0]]} ;
full_system.Vs{2} = {[[0,0];[0,1]],[[0]],[[1]]} ;
full_system.beta = beta ;

% information about the different blocks
full_system.block_coupling = struct() ;
full_system.block_coupling.E_blocks = [E_1,E_2,E_3] ;
full_system.block_coupling.coupled_blocks = [[1,2];[1,3];[2,3]] ;
full_system.block_coupling.coupling_matrices = {[Gamma_12;0],[0;Gamma_13],[Gamma_23]} ;
full_system.block_coupling.coupling_baths = {} ;
full_system.block_coupling.coupling_baths{1} = ...
    {struct("spectral_density","debye","omega_D",omega_12,"lambda_D",lambda_12,...
    "n_modes",n_modes)} ;
full_system.block_coupling.coupling_baths{2} = ...
    {struct("spectral_density","debye","omega_D",omega_13,"lambda_D",lambda_13,...
    "n_modes",n_modes)} ;
full_system.block_coupling.coupling_baths{3} = ...
    {struct("spectral_density","debye","omega_D",omega_23,"lambda_D",lambda_23,...
    "n_modes",n_modes)} ;
full_system.block_coupling.n_ts = [n_t,n_t,n_t] ;
full_system.block_coupling.t_maxs = [t_max,t_max,t_max] ;
% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition
heom_dynamics.rho_0_sys = {[[0,0];[0,1]],[[0]],[[0]]} ;

% set up observable arrays
heom_dynamics.observables = struct() ;
heom_dybamics.observables.block = {} ;
heom_dynamics.observables.block{1} = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;
heom_dynamics.observables.block{2} = {[[1]]} ;
heom_dynamics.observables.block{3} = {[[1]]} ;

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
% heom_dynamics.heom_truncation.truncation_method = "lambda weighted cut-off" ;
% heom_dynamics.heom_truncation.L_cut = L_cut ;
% heom_dynamics.heom_truncation.p = p ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;

% information on how to treat the block coupling
heom_dynamics.blocking_coupling = struct() ;
heom_dynamics.block_coupling.method = "truncated NZ" ;
heom_dynamics.block_coupling.Gamma_cut_trunc = Gamma_cut_trunc ;

% details of the strongly coupled bath
% AB_coupling_info = struct() ;
% AB_coupling_info.baths = {struct("spectral_density","debye","omega_D",omega_AB,"lambda_D",lambda_AB)} ;
% AB_coupling_info.beta = beta ;
% AB_coupling_info.coupling_matrix = [0;Gamma_AB] ; % Gamma_a,b matrix
% AB_coupling_info.Delta_E_AB = Delta_E_AB ;
% AB_coupling_info.t_max = t_max ;
% AB_coupling_info.n_t = n_t ;
% AB_coupling_info.n_modes = n_modes ;
% AB_coupling_info.method = "simplified" ;
% AB_coupling_info.Delta_E_AB = Delta_E_AB - (delta_epsilon/2);
% AB_coupling_info.method = "include H_sys" ;
% AB_coupling_info.method = "include H_sys NZ" ;
% AB_coupling_info.method = "full NZ" ;
% AB_coupling_info.method = "first-order phonon NZ" ;
% AB_coupling_info.method = "second-order phonon NZ" ;
% AB_coupling_info.method = "first-order phonon NZ 2" ;
% AB_coupling_info.method = "second-order phonon NZ 2" ;
% AB_coupling_info.method = "truncated NZ" ;
% AB_coupling_info.Gamma_cut = Gamma_cut_trunc ;


% run the dynamics
% [O_t_AB,t_AB,junk] = runHEOMTC2ABDynamics(full_system,heom_dynamics,AB_coupling_info) ;
% O_t_AB_full = O_t_AB ;
% skip = 1 ;
% t_AB_full = t_AB ;
% O_t_AB = O_t_AB(:,1:skip:end) ;
% t_AB = t_AB(1:skip:end) ;

% run the dynamics
[O_t_SCPT,t_SCPT,L,junk] = runHEOMSCPTDynamics(full_system,heom_dynamics) ;