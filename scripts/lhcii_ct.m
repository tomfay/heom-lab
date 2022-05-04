% This script tests the TCL2 approach for integrating out a strongly
% coupled bath from an exciton dimer coupled to a CT state

% Here the |A,k> states are |D*A> and |DA*> and the |B> state is |CT>.
% The |A,k> states are each coupled weakly to a debye bath.
% The |B> state is coupled strongly to a debye bath.

% set up system parameters
% H_s,A parameters
J_LE = [[0 36 -5 -3 1 -2 -3 3 4 -5 20 2 -8 2];
        [0 0 15 6 0 5 6 -6 -24 -5 1 8 -2 0 ];
        [0 0 0 -1 0 -4 6 4 72 7 -1 1 1 -5 ];
        [0 0 0 0 4 71 24 -4 -2 0 -3 3 2 -3 ];
        [0 0 0 0 0 9 -4 -4 0 1 1 -2 -1 0 ];
        [0 0 0 0 0 0 16 -5 2 0 -2 2 2 -2 ];
        [0 0 0 0 0 0 0 -4 -5 1 -2 3 3 -3 ];
        [0 0 0 0 0 0 0 0 24 43 5 -1 -2 1 ];
        [0 0 0 0 0 0 0 0 0 -2 4 -1 -2 2 ];
        [0 0 0 0 0 0 0 0 0 0 -26 13 6 -1 ];
        [0 0 0 0 0 0 0 0 0 0 0 99 -3 1 ];
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0 ];
        [0 0 0 0 0 0 0 0 0 0 0 0 0 -36]; 
        [0 0 0 0 0 0 0 0 0 0 0 0 0 0]] ;
J_LE = J_LE + transpose(J_LE) ;

E_LE = [15415
14850
14860
14920
15555
15395
15305
15175
15635
14780
14930
14960
14870
14980] ;
n_r = 12 ;
E_LE = E_LE - (E_LE(n_r)) ;
n_LE = length(E_LE) ;
H_sys_LE = J_LE + diag(E_LE) ;
[psi_exciton,E_exciton] = eig(H_sys_LE,'vector') ;

% set up explicit bath parameters
lambda_D = 220.0 ;
omega_D = 353.6777 ;
beta = 1.0/208.50907518 ;
lambda_AB = 5405.0 ;
omega_AB = 353.6777 ;
Gamma_AB = 0*240 ;
Delta_E_AB = 82.0 ;

% dynamics information
dt = 1e-3 ;
n_steps = 10000 ;
krylov_dim = 9 ;
krylov_tol = 1e-12 ;
Gamma_cut = 1200.0 ;
Gamma_cut_trunc = 00.00 ;
p = 1 ;
L_cut = 5.0 ;

% parameters for evaluating to AB correlation fction
t_max = sqrt((beta/lambda_AB)*log(1/1e-10)) ;
n_t = 4000 ;
n_modes = 512 ; % number of modes used to discretise the spectral density

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys_A = H_sys_LE ;
full_system.H_sys_B = [[0]] ;
% baths is a cell array of structs describign each bath
full_system.baths = {} ;
O_A = {} ;
for n = 1:n_LE
    V_A = sparse([n],[n],[1],n_LE,n_LE) ;
    O_A = [O_A,{V_A}] ;
    V_B = [[0]] ;
    if n == n_r 
        V_B = [[1]] ;
    end
    full_system.baths = [full_system.baths,...
    {struct("V_A",V_A,"V_B",V_B,...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)}] ;
end

% full_system.baths = [full_system.baths,...
%     {struct("V_A",[[1,0.5];[0.5,0]],"V_B",[[0]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",0.25*lambda_D)}] ;
% full_system.baths = {struct("V_A",[[1,0];[0,-1]],"V_B",[[0]],...
%     "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D)} ;
full_system.beta = beta ;

% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition
n_init = 9 ;
n_init_ex = 8 ;
heom_dynamics.rho_0_sys_A = zeros([n_LE,n_LE]) ;
heom_dynamics.rho_0_sys_A(n_init,n_init) = 1 ;
heom_dynamics.rho_0_sys_A = psi_exciton(:,n_init_ex)* psi_exciton(:,n_init_ex)' ;
heom_dynamics.rho_0_sys_B = [[0]] ;

% set up observable arrays
heom_dynamics.observables = struct() ;
coh = sparse([11],[12],[1],n_LE,n_LE) ;
O_A = [O_A,{coh}] ;
heom_dynamics.observables.system_A = O_A ;
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
% heom_dynamics.heom_truncation.truncation_method = "lambda weighted cut-off" ;
% heom_dynamics.heom_truncation.L_cut = L_cut ;
% heom_dynamics.heom_truncation.p = p ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;

% details of the strongly coupled bath
AB_coupling_info = struct() ;
AB_coupling_info.baths = {struct("spectral_density","debye","omega_D",omega_AB,"lambda_D",lambda_AB)} ;
AB_coupling_info.beta = beta ;
AB_coupling_info.coupling_matrix = sparse([n_r],[1],[Gamma_AB],n_LE,1) ; % Gamma_a,b matrix
AB_coupling_info.Delta_E_AB = Delta_E_AB ;
AB_coupling_info.t_max = t_max ;
AB_coupling_info.n_t = n_t ;
AB_coupling_info.n_modes = n_modes ;
% AB_coupling_info.method = "simplified" ;
% AB_coupling_info.method = "include H_sys" ;
% AB_coupling_info.method = "include H_sys NZ" ;
% AB_coupling_info.method = "full NZ" ;
% AB_coupling_info.method = "first-order phonon NZ" ;
% AB_coupling_info.method = "second-order phonon NZ" ;
% AB_coupling_info.method = "first-order phonon NZ 2" ;
% AB_coupling_info.method = "second-order phonon NZ 2" ;
AB_coupling_info.method = "truncated NZ" ;
AB_coupling_info.Gamma_cut = Gamma_cut_trunc ;


% run the dynamics
[O_t_AB,t_AB,junk] = runHEOMTC2ABDynamics(full_system,heom_dynamics,AB_coupling_info) ;
O_t_AB_full = O_t_AB ;
skip = 1 ;
t_AB_full = t_AB ;
O_t_AB = O_t_AB(:,1:skip:end) ;
t_AB = t_AB(1:skip:end) ;