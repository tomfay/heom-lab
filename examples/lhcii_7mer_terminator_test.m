% This script tests the TCL2 approach for integrating out a strongly
% coupled bath from an exciton dimer coupled to a CT state

% Here the |A,k> states are |D*A> and |DA*> and the |B> state is |CT>.
% The |A,k> states are each coupled weakly to a debye bath.
% The |B> state is coupled strongly to a debye bath.

% set up system parameters
% H_s,A parameters
H_sys_LE = [15889	49.64	-5.89	-2.51	0.77	-1.87	-2.49	2.78	3.79	-5.95	24.89	9.13	-10.79	3.59
49.64	15157	38.11	6.42	-0.71	5.6	7.13	-5.84	-19.25	-11.39	9.69	15.83	-4.96	0.69
-5.89	38.11	15287	-3.28	1.13	-8.89	1.23	6.72	96.66	12.97	-2.7	-0.76	2.68	-6.7
-2.51	6.42	-3.28	15460	3.35	104.56	35.93	-2.76	-7.28	-4.18	-3.8	4.67	2.12	-3.42
0.77	-0.71	1.13	3.35	15679	29.71	-4.47	-5.13	-0.77	1.61	1.33	-2.85	-1.4	0.37
-1.87	5.6	-8.89	104.56	29.71	15850	59.38	-4.99	-0.16	-3.28	-2.52	3.1	1.47	-2.16
-2.49	7.13	1.23	35.93	-4.47	59.38	15714	-4.43	-11.99	-0.14	-2.78	3.07	2.2	-3.25
2.78	-5.84	6.72	-2.76	-5.13	-4.99	-4.43	15761	36.07	61.97	4.35	-1.08	-2.01	1.3
3.79	-19.25	96.66	-7.28	-0.77	-0.16	-11.99	36.07	15721	3.86	4.3	-2.57	-2.92	2.33
-5.95	-11.39	12.97	-4.18	1.61	-3.28	-0.14	61.97	3.86	15073	-24.96	23.1	7.21	-1.55
24.89	9.69	-2.7	-3.8	1.33	-2.52	-2.78	4.35	4.3	-24.96	15112	126.92	-6.15	4.55
9.13	15.83	-0.76	4.67	-2.85	3.1	3.07	-1.08	-2.57	23.1	126.92	15094	-0.47	-0.18
-10.79	-4.96	2.68	2.12	-1.4	1.47	2.2	-2.01	-2.92	7.21	-6.15	-0.47	15174	-50.22
3.59	0.69	-6.7	-3.42	0.37	-2.16	-3.25	1.3	2.33	-1.55	4.55	-0.18	-50.22	15260] ;
inds_7mer = [2,3,8,9,10,11,12] ;
% inds_7mer = [8,9,10,11,12] ;
% inds_7mer = [10,11,12] ;
H_sys_LE = H_sys_LE(inds_7mer,inds_7mer) ;

E_LE = diag(H_sys_LE) ;

chlb_inds = find(E_LE>15000) ;
chla_inds = find(E_LE<=15000) ;
n_chla = numel(chla_inds) ;
n_LE = length(E_LE) ;
H_sys_LE = H_sys_LE -  eye(n_LE)*mean(E_LE) ;
[psi_exciton,E_exciton] = eig(H_sys_LE,'vector') ;

% set up explicit bath parameters
lambda_D_chla = 220.0 ;
omega_D_chla = 353.6777 ;
lambda_D_chlb = 220.0 ;
omega_D_chlb = 353.6777 ;
% lambda_D_chla = 37.0 ;
% omega_D_chla = 30.0 ;
% lambda_D_chlb = 48.0 ;
% omega_D_chlb = 30.0 ;
beta = 1.0/208.50907518 ;
beta = (300/150)*1.0/208.50907518 ; 


% dynamics information
dt = 1e-3 ;
n_steps = 4000 ;
krylov_dim = 16 ;
krylov_tol = 1e-10 ;
% Gamma_cut =  1* omega_D_chla ;
Gamma_cut = (2*pi/beta)*2.01 ;
L_max = 3 ;
M_max = 0 ;


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
% n_init_ex = 7 ;
heom_dynamics.rho_0_sys = psi_exciton(:,end)* psi_exciton(:,end)' ;
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
% heom_dynamics.heom_truncation.truncation_method = "depth cut-off" ;
% heom_dynamics.heom_truncation.M_max = M_max ;
% heom_dynamics.heom_truncation.L_max = L_max ;
% heom_dynamics.heom_truncation.heom_termination = "markovian" ;
heom_dynamics.heom_truncation.heom_termination = "low temp correction" ;
heom_dynamics.heom_truncation.heom_termination = "low temp correction RF2" ;
heom_dynamics.heom_truncation.heom_termination = "NZ2" ;
heom_dynamics.heom_truncation.heom_termination = "RF2" ;
% heom_dynamics.heom_truncation.diagonal_only_term = true ;
heom_dynamics.heom_truncation.termination_k_max = 50 ;
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