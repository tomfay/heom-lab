% set up system parameters
% H_s,A parameters
delta_epsilon = 1.0e0 ;
J = 1.0e0 ; 

% set up explicit bath parameters
lambda_D = 1.0e0 ;
omega_D = 1.75 ;
beta = 1.0 ;
lambda_AB = 5.0 ;
omega_AB = 1.75 ;
Gamma_AB = 0.5 ;
Delta_E_AB = 2.0 ;
eta = sqrt(lambda_AB/lambda_D) ;

% dynamics information
dt = 0.01e0 ;
n_steps = 10000 ;
krylov_dim = 9 ;
krylov_tol = 1e-12 ;
order_adapt_taylor = 4 ;
tol_adapt_taylor = 1e-2 ;
Gamma_cut = 5 ;
p = 1 ;
L_cut = 40 ;

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = [[delta_epsilon/2,J,0];
                       [J,-delta_epsilon/2,Gamma_AB];
                       [0,Gamma_AB,lambda_AB-Delta_E_AB]];
% baths is a cell array of structs describign each bath
full_system.baths = {} ;
full_system.baths = [full_system.baths,...
    {struct("V",[[1/sqrt(2),0,0];[0,-1/sqrt(2),0];[0,0,-1/sqrt(2)]],...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)}] ;
full_system.baths = [full_system.baths,...
    {struct("V",[[0,0,0];[0,0,0];[0,0,1]],...
    "spectral_density","debye","omega_D",omega_AB,"lambda_D",lambda_AB)}] ;
full_system.beta = beta ;

% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition
heom_dynamics.rho_0_sys = [[1,0,0];[0,0,0];[0,0,0]] ;


% set up observable arrays
heom_dynamics.observables = struct() ;
heom_dynamics.observables.system = {[[0,1,0];[1,0,0];[0,0,0]],[[0,-1.0i,0];[1.0i,0,0];[0,0,0]],[[1,0,0];[0,-1,0];[0,0,0]],[[1,0,0];[0,1,0];[0,0,0]],[[0,0,0];[0,0,0];[0,0,1]]} ;

% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct() ;
heom_dynamics.integrator.method = "adaptive SIA" ;
heom_dynamics.integrator.dt = dt ;
heom_dynamics.integrator.n_steps = n_steps ;
heom_dynamics.integrator.krylov_dim = krylov_dim ;
heom_dynamics.integrator.krylov_tol = krylov_tol ;
% heom_dynamics.integrator = struct() ;
% heom_dynamics.integrator.method = "adaptive taylor" ;
% heom_dynamics.integrator.order = order_adapt_taylor ;
% heom_dynamics.integrator.tol = tol_adapt_taylor ;
% heom_dynamics.integrator.t_max = n_steps*dt ;

% hierarchy trunction information
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
% heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "coupling weighted cut-off" ;
heom_dynamics.heom_truncation.truncation_method = "lambda weighted cut-off" ;
heom_dynamics.heom_truncation.L_cut = L_cut ;
heom_dynamics.heom_truncation.p = p ;
heom_dynamics.heom_truncation.heom_termination = "markovian" ;


% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
O_t_full = O_t ;
t_full = t ;
skip = 1 ;
O_t = O_t(:,1:skip:end) ;
t = t(1:skip:end) ;