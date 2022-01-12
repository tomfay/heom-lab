% test the HEOM stuff with a 2 level system coupled to a bath
Delta = 0.0 ;
beta = 1.0 ;
gamma_B = [0.2,0.2] ;
lambda_B = [1.0,1.0] ;
Omega_B = [2.0,2.0]  ;
% gamma_B = [0.2] ;
% Omega_B = [2.0] ;
% lambda_B = [1.0]s ;
epsilon = 10.0 ;
k_B_T = 1/beta ;

% define the system Hamiltonian and the system-bath coupling
d_hilb = 2 ;
d_liou = d_hilb^2 ;
H_sys = sparse([0,Delta;Delta,epsilon+lambda_B(1)]) ;
V = {sparse([0,0;0,1]),sparse([0,1;1,0])} ;
% V = {sparse([0,0;0,1])} ;


% construct the system liouvillian
id_sys = speye(2) ;
L_sys = -1.0i * (kron(H_sys,id_sys)-kron(id_sys,transpose(H_sys))) ;

% initial system state
% rho_0_sys = [1;0;0;0] ;
rho_0_sys = [0;0;0;1] ;
% rho_0_sys = [0;1;1;0] ;
rho_0_sys = convertToLiouvilleVector([[0,-1i];[1i,0]]) ;

% build the HEOM generator, setting the truncation level
% N_trunc = 100 ; 
% L_HEOM = constructHiT1BathHEOMGenerator(L_sys,V,beta,gamma,lambda,N_trunc) ;
Gamma_cut = 8.0 * Omega_B(1) ;
fprintf('Creating HEOM generator... ') ;
[L_heom,ado_indices] = constructScaledHEOMGeneratorFreqTruncUBO(H_sys,V,gamma_B,lambda_B,Omega_B,beta,Gamma_cut) ;
n_ados = size(ado_indices,1) ;
d_heom = d_liou*n_ados ;
fprintf('Done.\nn_ados = %i\n',n_ados) ;
rho_0_heom = zeros([(n_ados)*d_liou,1]) ;
rho_0_heom(1:d_liou) = rho_0_sys ;

% observables
O = sparse([],[],[],4,(n_ados)*d_liou) ;
O(1,1:d_liou) = [1 , 0, 0, 0 ] ;
O(2,1:d_liou) = [0 , 0, 0, 1 ] ;
O(3,1:d_liou) = [0, 1, 1, 0] ;
O(4,1:d_liou) = [0, -1, 1, 0] ;

% calculate the marcus rate
% k_marcus = 2*pi*(Delta^2)/sqrt(4*pi*lambda_B*k_B_T)*exp(-(lambda_B + epsilon)^2 / (4*lambda_B * k_B_T)) ;
% t_tot = 100/k_marcus ;
t_tot = 100 ;


% dynamics with SIA
dim_krylov = 16 ;
tol = 0.25e-10 ;
% dt = 0.5e-1 ;
n_steps = 40000 ;
dt = t_tot / n_steps ;
fprintf('Starting dynamics... ');
[O_t,t] = runDynamicsSIADensityOperator(rho_0_heom,L_heom,n_steps,dt,O,dim_krylov,tol) ;
fprintf('Done.\n');

% t_heom = t ;
% % t_heom = 0:0.01:10 ;
% dt_heom = t_heom(2)- t_heom(1) ;
% omega_fft = (2 * pi/(length(t_heom)*dt_heom)) * (0:(length(t_heom)-1)) ;
% omega_fft = -(omega_fft - omega_fft(end)/2) ;
% ft_cutoff_func = (0.5+0.5*cos(t*pi/t(end))) ;
% chi = dt_heom*fftshift(fft(O_t(3,:).*ft_cutoff_func)) ;
% chi_imag = imag(chi) ;

n_t = length(t) ;
n_omega = 4*n_t ;
ft_cutoff_func = (0.5+0.5*cos(t*pi/t(end))) ;
S_t = [(O_t(3,:)).*ft_cutoff_func,zeros(1,n_omega-n_t)];
t_pad = dt*(0:(n_omega-1)) ;


omega_fft = (2 * pi/(length(t_pad)*dt)) * (0:(length(t_pad)-1)) ;
omega_fft = -(omega_fft - omega_fft(end)/2) ;
chi = dt*fftshift(fft(S_t)) ;
chi_imag = imag(chi) ;

% Omega_test = 2 ; tau_test = 2 ;
% f_test = exp(1.0i*Omega_test*t_heom - (t_heom/tau_test).^2) ;
% ft_f_test = 0.5*sqrt(pi)*tau_test*exp(-0.25*(tau_test^2)*(omega_fft+Omega_test).^2)...
%     + 1.0i*tau_test * dawson(0.5*tau_test*(omega_fft+Omega_test));
% % .*(1+1.0i*erfi(0.5*tau_test*(omega_fft+Omega_test)));
% fft_f_test = dt_heom * fftshift(fft(f_test)) ;

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

