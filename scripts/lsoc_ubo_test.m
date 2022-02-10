% A test script for the spin boson model
% In this example the dynamics for a spin boson model are calculated
% In this example the spin is coupled to two baths a debye bath and an
% underdamped brownian oscillator bath

% Parameters for the problem
% system hamiltonian parameters
epsilon = 0.75 ;
% bath parameters
beta = 1.0 ;
% debye bath parameters
lambda_D_t = 0.6 ;
omega_D_t = 0.2 ;
lambda_c = 1.0 ;
Omega_c = 1.0 ;
gamma_c = 0.2 ;
Delta = 1.0 ;
Lambda = 0.1 ;
phi_0 = pi*0.5 ;
Delta_0 = Delta * cos(phi_0) ;
Delta_1 = -Delta * sin(phi_0) ;
Lambda_0 = Lambda * sin(phi_0) ;
Lambda_1 = Lambda * cos(phi_0) ;
eta = 0.0 ; % determines how close to eq initial density is

% dynamics information
dt = 1e-2 ;
n_steps = 10000 ;
krylov_dim = 8 ;
krylov_tol = 1e-8 ;
Gamma_cut = 8.0 ;

% matrices of system observable operators to be returned, sigma_x, sigma_y
% sigma_z, and 1
O_sys = {[[0,1];[1,0]],[[0,-1.0i];[1.0i,0]],[[1,0];[0,-1]],eye(2)} ;

% initial state of the system
rho_0_sys = [[1,0];[0,0]] ;

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon/2,Delta_0+1.0i*Lambda_0];
                     [Delta_0-1.0i*Lambda_0,-epsilon/2]];
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",[[1,0];[0,-1]],...
    "spectral_density","debye","omega_D",omega_D_t,"lambda_D",lambda_D_t)} ;
full_system.baths = [full_system.baths,struct("V",[[1*(1-eta),(Delta_1 + 1.0i * Lambda_1)];[(Delta_1 - 1.0i * Lambda_1),-1*(1+eta)]],...
    "spectral_density","UBO","Omega",Omega_c,"lambda",lambda_c,...
    "gamma",gamma_c)] ;
full_system.beta = beta ;

% a struct that contains information about the HEOM dynamics
heom_dynamics = struct ;
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
heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
heom_dynamics.heom_truncation.heom_termination = "markovian" ;

% what system observables should be returned
heom_dynamics.observables = struct ;
heom_dynamics.observables.system = O_sys ;

% set the initial condition
heom_dynamics.rho_0_sys = rho_0_sys ;

% run the dynamics
[O_t_up,t] = runHEOMDynamics(full_system,heom_dynamics) ;

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[epsilon/2,Delta_0-1.0i*Lambda_0];
                     [Delta_0+1.0i*Lambda_0,-epsilon/2]];
% baths is a cell array of structs describign each bath
full_system.baths = {struct("V",[[1,0];[0,-1]],...
    "spectral_density","debye","omega_D",omega_D_t,"lambda_D",lambda_D_t)} ;
full_system.baths = [full_system.baths,struct("V",[[1*(1-eta),(Delta_1 - 1.0i * Lambda_1)];[(Delta_1 + 1.0i * Lambda_1),-1*(1+eta)]],...
    "spectral_density","UBO","Omega",Omega_c,"lambda",lambda_c,...
    "gamma",gamma_c)] ;
full_system.beta = beta ;

% run the dynamics
[O_t_down,t] = runHEOMDynamics(full_system,heom_dynamics) ;