% This script tests the TCL2 approach for integrating out a strongly
% coupled bath from an exciton dimer coupled to a CT state

% Here the |A,k> states are |D*A> and |DA*> and the |B> state is |CT>.
% The |A,k> states are each coupled weakly to a debye bath.
% The |B> state is coupled strongly to a debye bath.

% set up system parameters
% H_s,A parameters
J_LE = [0 3.811000e+01 -1.139000e+01 9.690000e+00 1.583000e+01 -5.840000e+00 -1.925000e+01 -4.960000e+00 6.900000e-01 6.420000e+00 -7.100000e-01 5.600000e+00 7.130000e+00 4.964000e+01 
3.811000e+01 0 1.297000e+01 -2.700000e+00 -7.600000e-01 6.720000e+00 9.666000e+01 2.680000e+00 -6.700000e+00 -3.280000e+00 1.130000e+00 -8.890000e+00 1.230000e+00 -5.890000e+00 
-1.139000e+01 1.297000e+01 0 -2.496000e+01 2.310000e+01 6.197000e+01 3.860000e+00 7.210000e+00 -1.550000e+00 -4.180000e+00 1.610000e+00 -3.280000e+00 -1.400000e-01 -5.950000e+00 
9.690000e+00 -2.700000e+00 -2.496000e+01 0 1.269200e+02 4.350000e+00 4.300000e+00 -6.150000e+00 4.550000e+00 -3.800000e+00 1.330000e+00 -2.520000e+00 -2.780000e+00 2.489000e+01 
1.583000e+01 -7.600000e-01 2.310000e+01 1.269200e+02 0 -1.080000e+00 -2.570000e+00 -4.700000e-01 -1.800000e-01 4.670000e+00 -2.850000e+00 3.100000e+00 3.070000e+00 9.130000e+00 
-5.840000e+00 6.720000e+00 6.197000e+01 4.350000e+00 -1.080000e+00 0 3.607000e+01 -2.010000e+00 1.300000e+00 -2.760000e+00 -5.130000e+00 -4.990000e+00 -4.430000e+00 2.780000e+00 
-1.925000e+01 9.666000e+01 3.860000e+00 4.300000e+00 -2.570000e+00 3.607000e+01 0 -2.920000e+00 2.330000e+00 -7.280000e+00 -7.700000e-01 -1.600000e-01 -1.199000e+01 3.790000e+00 
-4.960000e+00 2.680000e+00 7.210000e+00 -6.150000e+00 -4.700000e-01 -2.010000e+00 -2.920000e+00 0 -5.022000e+01 2.120000e+00 -1.400000e+00 1.470000e+00 2.200000e+00 -1.079000e+01 
6.900000e-01 -6.700000e+00 -1.550000e+00 4.550000e+00 -1.800000e-01 1.300000e+00 2.330000e+00 -5.022000e+01 0 -3.420000e+00 3.700000e-01 -2.160000e+00 -3.250000e+00 3.590000e+00 
6.420000e+00 -3.280000e+00 -4.180000e+00 -3.800000e+00 4.670000e+00 -2.760000e+00 -7.280000e+00 2.120000e+00 -3.420000e+00 0 3.350000e+00 1.045600e+02 3.593000e+01 -2.510000e+00 
-7.100000e-01 1.130000e+00 1.610000e+00 1.330000e+00 -2.850000e+00 -5.130000e+00 -7.700000e-01 -1.400000e+00 3.700000e-01 3.350000e+00 0 2.971000e+01 -4.470000e+00 7.700000e-01 
5.600000e+00 -8.890000e+00 -3.280000e+00 -2.520000e+00 3.100000e+00 -4.990000e+00 -1.600000e-01 1.470000e+00 -2.160000e+00 1.045600e+02 2.971000e+01 0 5.938000e+01 -1.870000e+00 
7.130000e+00 1.230000e+00 -1.400000e-01 -2.780000e+00 3.070000e+00 -4.430000e+00 -1.199000e+01 2.200000e+00 -3.250000e+00 3.593000e+01 -4.470000e+00 5.938000e+01 0 -2.490000e+00 
4.964000e+01 -5.890000e+00 -5.950000e+00 2.489000e+01 9.130000e+00 2.780000e+00 3.790000e+00 -1.079000e+01 3.590000e+00 -2.510000e+00 7.700000e-01 -1.870000e+00 -2.490000e+00 0] ;


E_LE = [15157
15287
15073
15112
15094
15761
15721
15174
15260
15460
15679
15850
15714
15889] ;
chlb_inds = find(E_LE>15500) ;
chla_inds = find(E_LE<=15500) ;
n_612a = 5 ;
n_612a = find(chla_inds==n_612a) ;
n_603a = 2 ;
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
lambda_D = 37.0 ;
omega_D = 30.0 ;
beta = 1.0/208.50907518 ;
lambda_lut1CT = 5405.0 ;
lambda_lut2CT = 5052.0 ;
omega_lut1CT = omega_D ;
omega_lut2CT = omega_D ;
Gamma_lut1CT = 240 ;
Gamma_lut2CT = 279 ;
E_lut1CT = E_LE(n_612a)-82 ;
E_lut2CT = E_LE(n_603a)+951 ;
E_gs = -E_LE_max+lambda_D ;

% dynamics information
dt = 1e-3 ;
n_steps = 50000 ;
krylov_dim = 16 ;
krylov_tol = 1e-8 ;
Gamma_cut = 2.1 * omega_D ;
Gamma_cut_trunc = 1.1 *omega_D ;
p = 1 ;
L_cut = 5.0 ;

% parameters for evaluating to AB correlation fction
t_max = sqrt((beta/lambda_lut1CT)*log(1/1e-10)) ;
n_t = 1000 ;
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
% full_system.block_coupling.coupling_baths{1} = ...
%     {struct("spectral_density","high temperature","lambda",lambda_lut1CT)} ;
% full_system.block_coupling.coupling_baths{2} = ...
%     {struct("spectral_density","high temperature","lambda",lambda_lut2CT)} ;
% full_system.block_coupling.coupling_baths{3} = ...
%     {struct("spectral_density","high temperature","lambda",lambda_lut1CT)} ;
% full_system.block_coupling.coupling_baths{4} = ...
%     {struct("spectral_density","high temperature","lambda",lambda_lut2CT)} ;

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
heom_dynamics.integrator.method = "adaptive SIA" ;
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
c_0 = 2.99792458e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;