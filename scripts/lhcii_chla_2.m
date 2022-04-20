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
14980 ] ;
chlb_inds = find(E_LE>15000) ;
chla_inds = find(E_LE<=15000) ;
n_612a = 12 ;
n_612a = find(chla_inds==n_612a) ;
n_603a = 3 ;
n_603a = find(chla_inds==n_603a) ;
% only include Chlorophyll A
E_LE = E_LE(chla_inds) ;
E_LE_max = max(E_LE) ;
E_LE = E_LE - max(E_LE) ;
J_LE = J_LE(chla_inds(1:(end)),chla_inds(1:(end))) ;
n_LE = length(E_LE) ;
H_sys_LE = diag(E_LE) ;
H_sys_LE(1:(n_LE),1:(n_LE)) = H_sys_LE(1:(n_LE),1:(n_LE)) + J_LE ;
[psi_exciton,E_exciton] = eig(H_sys_LE,'vector') ;

% set up explicit bath parameters
lambda_D = 220.0 ;
omega_D = 353.6777 ;
beta = 1.0/208.50907518 ;
lambda_lut1CT = 5405.0 ;
lambda_lut2CT = 5052.0 ;
omega_lut1CT = 353.6777*0.25 ;
omega_lut2CT = 353.6777*0.25 ;
Gamma_lut1CT = 240 ;
Gamma_lut2CT = 279 ;
E_lut1CT = E_LE(n_612a)-82 ;
E_lut2CT = E_LE(n_603a)+951 ;
E_gs = -E_LE_max+lambda_D ;

% dynamics information
dt = 1e-3 ;
n_steps = 50000 ;
krylov_dim = 9 ;
krylov_tol = 1e-12 ;
Gamma_cut = 1200.0 ;
Gamma_cut_trunc = 00.00 ;
p = 1 ;
L_cut = 5.0 ;

% parameters for evaluating to AB correlation fction
t_max = sqrt((beta/lambda_lut1CT)*log(1/1e-10)) ;
n_t = 400 ;
n_modes = 512 ; % number of modes used to discretise the spectral density

% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = {} ;
full_system.H_sys{1} = H_sys_LE ;
full_system.H_sys{2} = [[0]] ;
full_system.H_sys{3} = [[0]] ;
full_system.H_sys{4} = [[0]] ;
% baths is a cell array of structs describign each bath
full_system.baths = cell([n_LE,1]) ;
full_system.Vs = cell([n_LE,1]) ;
O_LE = {} ;
for n = 1:(n_LE)
    V_LE = sparse([n],[n],[1],n_LE,n_LE) ;
    O_LE = [O_LE,{V_LE}] ;
    V_lut1CT = [[0]] ;
    V_lut2CT = [[0]] ;
    V_gs = [[0]] ;
    if n == n_612a  
        V_lut1CT = [[1]] ;
    end
    if n == n_603a  
        V_lut1CT = [[1]] ;
    end
    full_system.baths{n} = struct(...
    "spectral_density","debye","omega_D",omega_D,"lambda_D",lambda_D) ;
    full_system.Vs{n} = {V_LE,V_lut1CT,V_lut2CT,V_gs} ;
end
% baths is a cell array of structs describign each bath
full_system.beta = beta ;

% information about the different blocks
full_system.block_coupling = struct() ;
full_system.block_coupling.E_blocks = [0,E_lut1CT,E_lut2CT,E_gs] ;
full_system.block_coupling.coupled_blocks = [[1,2];[1,3];[2,4];[3,4]] ;
full_system.block_coupling.coupling_matrices = {zeros([n_LE,1]),zeros([n_LE,1]),[[Gamma_lut1CT]],[[Gamma_lut2CT]]} ;
full_system.block_coupling.coupling_matrices{1}(n_612a) = Gamma_lut1CT ;
full_system.block_coupling.coupling_matrices{2}(n_603a) = Gamma_lut2CT ;
full_system.block_coupling.coupling_baths = {} ;
full_system.block_coupling.coupling_baths{1} = ...
    {struct("spectral_density","debye","omega_D",omega_lut1CT,"lambda_D",lambda_lut1CT,...
    "n_modes",n_modes)} ;
full_system.block_coupling.coupling_baths{2} = ...
    {struct("spectral_density","debye","omega_D",omega_lut2CT,"lambda_D",lambda_lut2CT,...
    "n_modes",n_modes)} ;
full_system.block_coupling.coupling_baths{3} = ...
    {struct("spectral_density","debye","omega_D",omega_lut1CT,"lambda_D",lambda_lut1CT,...
    "n_modes",n_modes)} ;
full_system.block_coupling.coupling_baths{4} = ...
    {struct("spectral_density","debye","omega_D",omega_lut2CT,"lambda_D",lambda_lut2CT,...
    "n_modes",n_modes)} ;
full_system.block_coupling.n_ts = [n_t,n_t,n_t,n_t] ;
full_system.block_coupling.t_maxs = [t_max,t_max,t_max,t_max] ;

% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition
n_init_ex = 8 ;
heom_dynamics.rho_0_sys = {psi_exciton(:,n_init_ex)* psi_exciton(:,n_init_ex)',[[0]],[[0]],[[0]]} ;
heom_dynamics.rho_0_sys = {(1/n_LE)*eye(n_LE),[[0]],[[0]],[[0]]} ;
% set up observable arrays
heom_dynamics.observables = struct() ;
% coh = sparse([11],[12],[1],n_LE,n_LE) ;
% O_LE = [O_LE,{coh}] ;
heom_dynamics.observables.block = {} ;
heom_dynamics.observables.block{1} = O_LE ;
heom_dynamics.observables.block{2} = {[[1]]} ;
heom_dynamics.observables.block{3} = {[[1]]} ;
heom_dynamics.observables.block{4} = {[[1]]} ;

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
% information on how to treat the block coupling
heom_dynamics.blocking_coupling = struct() ;
heom_dynamics.block_coupling.method = "truncated NZ" ;
heom_dynamics.block_coupling.Gamma_cut_trunc = Gamma_cut_trunc ;


% run the dynamics
[O_t,t,L,junk] = runHEOMSCPTDynamics(full_system,heom_dynamics) ;
c_0 = 3e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;