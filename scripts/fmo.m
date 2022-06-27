% This script tests the TCL2 approach for integrating out a strongly
% coupled bath from an exciton dimer coupled to a CT state

% Here the |A,k> states are |D*A> and |DA*> and the |B> state is |CT>.
% The |A,k> states are each coupled weakly to a debye bath.
% The |B> state is coupled strongly to a debye bath.

% set up system parameters
% H_s,A parameters
H_sys_LE = [410	-87.7	5.5	-5.9	6.7	-13.7	-9.9
-87.7	530	30.8	8.2	0.7	11.8	4.3
5.5	30.8	210	-53.5	-2.2	-9.6	6.0
-5.9	8.2	-53.5	320	-70.7	-17.0	-63.3
6.7	0.7	-2.2	-70.7	480	81.1	-1.3
-13.7	11.8	-9.6	-17.0	81.1	630	39.7
-9.9	4.3	6.0	-63.3	-1.3	39.7	440] ;
% inds_7mer = [2,3,8,9,10,11,12] ;
% inds_7mer = [8,9,10,11,12] ;
% inds_7mer = [10,11,12] ;
% H_sys_LE = H_sys_LE(inds_7mer,inds_7mer) ;

E_LE = diag(H_sys_LE) ;

chlb_inds = find(E_LE>15000) ;
chla_inds = find(E_LE<=15000) ;
n_chla = numel(chla_inds) ;
n_LE = length(E_LE) ;
H_sys_LE = H_sys_LE -  eye(n_LE)*mean(E_LE) ;
[psi_exciton,E_exciton] = eig(H_sys_LE,'vector') ;

% set up explicit bath parameters
lambda_D_chla = 35 ;
omega_D_chla = 106.1 ;
lambda_D_chlb = 35 ;
omega_D_chlb = 106.1 ;
% omega_D_chla = 106.1*(166/50) ;
% omega_D_chlb = 106.1*(166/50) ;
% lambda_D_chla = 37.0 ;
% omega_D_chla = 30.0 ;
% lambda_D_chlb = 48.0 ;
% omega_D_chlb = 30.0 ;
beta = 1.0/208.50907518 ;
% beta = (300/77)*1.0/208.50907518 ; 


% dynamics information
dt = 1.0e-3 ;
n_steps = 200 ;
krylov_dim = 16 ;
krylov_tol = 1e-10 ;
% Gamma_cut =  3.01* omega_D_chla ;
Gamma_cut = (1*pi/beta)*2.01 ;
L_max = 3 ;
M_max = 1 ;


% the full_system object contains all information about the Hamiltonian of
% the full open quantum system
full_system = struct() ;
% H_sys contains the system Hamiltonian
full_system.H_sys = H_sys_LE ;
% baths is a cell array of structs describign each bath
full_system.beta = beta ;
% baths is a cell array of structs describign each bath
full_system.baths = cell([n_LE,1]) ;
O_LE = {} ;
for n = 1:(n_LE)
    V_LE = sparse([n],[n],[1],n_LE,n_LE) ;
    O_LE = [O_LE,{V_LE}] ;
    if ismember(n,chla_inds)
        full_system.baths{n} = struct("V",V_LE,...
            "spectral_density","debye","omega_D",omega_D_chla,"lambda_D",lambda_D_chla) ;
    else
        full_system.baths{n} = struct("V",V_LE,...
            "spectral_density","debye","omega_D",omega_D_chlb,"lambda_D",lambda_D_chlb) ;
    end
end
% baths is a cell array of structs describign each bath
full_system.beta = beta ;


% set up the dynamics struct
heom_dynamics = struct() ;
% set the initial condition

heom_dynamics.rho_0_sys = psi_exciton(:,end)* psi_exciton(:,end)' ;
n_init_site = 1 ;
heom_dynamics.rho_0_sys = zeros([n_LE,n_LE]) ;
heom_dynamics.rho_0_sys(n_init_site,n_init_site) = 1.0 ;

% heom_dynamics.rho_0_sys{1} = expm(-beta*H_sys_LE)/sum(diag(expm(-beta*H_sys_LE))) ;
% set up observable arrays
heom_dynamics.observables = struct() ;
% coh = sparse([11],[12],[1],n_LE,n_LE) ;
% O_LE = [O_LE,{coh}] ;
heom_dynamics.observables.system = O_LE ;

% integrator information, currently only the short iterative arnoldi is
% implemented
heom_dynamics.integrator = struct() ;
heom_dynamics.integrator.method = "adaptive SIA" ;
heom_dynamics.integrator.dt = dt ;
heom_dynamics.integrator.n_steps = n_steps ;
heom_dynamics.integrator.krylov_dim = krylov_dim ;
heom_dynamics.integrator.krylov_tol = krylov_tol ;

% hierarchy trunction information
heom_dynamics.heom_truncation = struct() ;
heom_dynamics.heom_truncation.truncation_method = "frequency cut-off" ;
heom_dynamics.heom_truncation.Gamma_cut = Gamma_cut ;
heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
heom_dynamics.heom_truncation.M_max = M_max ;
heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
heom_dynamics.heom_truncation.heom_termination = "low temp correction" ;
% heom_dynamics.heom_truncation.heom_termination = "low temp correction NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "NZ2" ;
% heom_dynamics.heom_truncation.heom_termination = "RF2" ;
heom_dynamics.heom_truncation.diagonal_only_term = true ;
heom_dynamics.heom_truncation.termination_k_max = 200 ;
% heom_dynamics.heom_truncation = struct() ;
% heom_dynamics.heom_truncation.truncation_method = "lambda weighted cut-off" ;
% heom_dynamics.heom_truncation.L_cut = L_cut ;
% heom_dynamics.heom_truncation.p = p ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;

% run the dynamics
[O_t,t] = runHEOMDynamics(full_system,heom_dynamics) ;
n_t_plot = min([1000,n_steps]) ;
skip = floor(n_steps/n_t_plot);
O_t = O_t(:,1:skip:end) ;
t = t(1:skip:end) ;
c_0 = 2.99792458e10 ;
t_ps = t / (2*pi*c_0 * 1e-12) ;