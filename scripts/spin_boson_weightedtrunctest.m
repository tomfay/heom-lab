% test the HEOM stuff with a 2 level system coupled to a bath
epsilon = 0.0 ;
Delta = 2.0 ;
beta = 1.0 ;
gamma_D = [1.0,1.0] ;
lambda_D = [5.0,0.5] ;
k_B_T = 1/beta ;

% define the system Hamiltonian and the system-bath coupling
d_hilb = 2 ;
d_liou = d_hilb^2 ;
H_sys = sparse([0,Delta;Delta,epsilon+lambda_D(1)]) ;
V = {sparse([0,0;0,1]),sparse([0,1;1,0])} ;

% construct the system liouvillian
id_sys = speye(2) ;
L_sys = -1.0i * (kron(H_sys,id_sys)-kron(id_sys,transpose(H_sys))) ;

% initial system state
rho_0_sys = [1;0;0;0] ;

% build the HEOM generator, setting the truncation level
% N_trunc = 100 ; 
% L_HEOM = constructHiT1BathHEOMGenerator(L_sys,V,beta,gamma,lambda,N_trunc) ;
Gamma_cut = 30.0 * gamma_D(1) ;
cutoff_param = 17.0 ;
fprintf('Creating HEOM generator... ') ;
% [L_heom,ado_indices] = constructScaledHEOMGeneratorFreqTrunc(H_sys,V,gamma_D,lambda_D,beta,Gamma_cut) ;
[L_heom,ado_indices] = constructScaledHEOMGeneratorWeightedCutoff(H_sys,V,gamma_D,lambda_D,beta,cutoff_param) ;

n_ados = size(ado_indices,1) ;
d_heom = d_liou*n_ados ;
fprintf('Done.\nn_ados = %i\n',n_ados) ;
rho_0_heom = zeros([(n_ados)*d_liou,1]) ;
rho_0_heom(1:d_liou) = rho_0_sys ;

% observables
O = sparse([],[],[],2,(n_ados)*d_liou) ;
O(1,1:d_liou) = [1 , 0, 0, 0 ] ;
O(2,1:d_liou) = [0 , 0, 0, 1 ] ;

% calculate the marcus rate
% k_marcus = 2*pi*(Delta^2)/sqrt(4*pi*lambda_D*k_B_T)*exp(-(lambda_D + epsilon)^2 / (4*lambda_D * k_B_T)) ;
% t_tot = 100/k_marcus ;
t_tot = 50 ;


% dynamics with SIA
dim_krylov = 16 ;
tol = 1e-13 ;
% dt = 0.5e-1 ;
n_steps = 100000 ;
dt = t_tot / n_steps ;
fprintf('Starting dynamics... ');
[O_t,t] = runDynamicsSIADensityOperator(rho_0_heom,L_heom,n_steps,dt,O,dim_krylov,tol) ;
fprintf('Done.\n');

% % calculate the kernel
% A_rho = O' ;
% A_rho_sp = sparse([],[],[],d_heom,2) ;
% A_rho_sp(1:d_liou,1:2) = O(1:2,1:d_liou)' ;
% A = O ;
% A_sp = sparse([],[],[],2,d_heom) ;
% A_sp(1:2,1:d_liou) = O(1:2,1:d_liou) ;
% % P = A_rho * A ;
% P = A_rho_sp * A_sp ;
% Q = speye(d_heom) - P ;
% QL = Q * L_heom ;
% ALQ  = A * L_heom *  Q ;
% n_A = size(A_rho,2) ;
% kappa = zeros(n_A^2,n_steps+1 ) ;
% for n = 1:n_A
%     rho_0 = QL*A_rho(:,n) ;
%     [kernels_t,t] = runDynamicsSIADensityOperator(rho_0,QL,n_steps,dt,ALQ,dim_krylov,tol) ;
%     kappa(((n-1)*n_A+1):(n*n_A),:) = kernels_t ;
% end
% PLP = A*L_heom*A_rho ;
% n_mem = floor((n_steps+1)) ;
% kappa_trunc = kappa(:,1:n_mem) ;
% K = reshape(kappa_trunc,[2,2*size(kappa_trunc,2)]) ;
% b = zeros([2,3]) ;
% y_0 = [1;0] ;
% y_t = integrateIntegroDiffVolterraEquation(y_0,dt,PLP,b,real(K),n_steps,'trapz') ;
% y_t_simp = integrateIntegroDiffVolterraEquation(y_0,dt,PLP,b,real(K),n_steps,'simps') ;
% y_t_simps_alt = integrateIntegroDiffVolterraEquation(y_0,dt,PLP,b,real(K),n_steps,'simps-alt') ;