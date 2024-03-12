
%bias
epsilon = 93.3658 ;
epsilon  = 20 ;
% epsilon = 9.2662e+01 ;
% bath parameters
beta = 1 ;
% debye bath parameterslogl
% debye bath parameters
lambda_D = (11.0291)*0 ;
lambda_D = 10.1459*0 ;
omega_D =  0.1831 ;
% omega_D = 0.5 ;
% lambda_D = 20 ;
lambda_D = 20 ; omega_D = 0.5 ;
% UBO bath parameters Omega > gamma/2
% Omega_UBO = [6.7600;0.25] ;
% gamma_UBO = [2;0.1]  ;
% lambda_UBO = [8.6780;10.1459*0] ; 
Omega_UBO = [1] ;
gamma_UBO = [16]  ;
lambda_UBO = [0] ; 
% lambda_UBO = [3.9633e+00] ;
% lambda_UBO = 8.8830e+00 ;
Omega_OBO = [] ;
gamma_OBO = []  ;
lambda_OBO = [] ; 
%combine UBO and OBO
Omega_UBO = [Omega_UBO;Omega_OBO] ;
gamma_UBO = [gamma_UBO;gamma_OBO] ;
lambda_UBO = [lambda_UBO;lambda_OBO] ;
% total reorg
Lambda = sum(lambda_D) + sum(lambda_UBO) ;

% omega = 0:0.1:10 ;
Omega_max = 250 ;
omega = [linspace(0,Omega_max,2000)] ;

j_D = sum((0.5.*lambda_D) .* omega_D ./ (omega_D.^2 + omega.^2) ,1) ;
% j_UBO = (0.5*lambda_UBO) * Omega_UBO^2 * gamma_UBO ./ ((omega.^2-Omega_UBO^2).^2 + gamma_UBO^2 * omega.^2) ;

j_UBO = sum((0.5.*lambda_UBO) .* Omega_UBO.^2 .* gamma_UBO ./ ( (omega.^2-Omega_UBO.^2).^2 + gamma_UBO.^2 .* omega.^2),1) ;
% j_UBO =  lambda_UBO * sqrt(1/(2*pi*gamma_UBO^2))*exp(-(omega-Omega_UBO).^2 / (2*gamma_UBO.^2)) ;

j_omega_eff = j_D + j_UBO ;




norm_rho = 1/trapz(omega,j_omega_eff) ;
j_omega_eff = j_omega_eff * norm_rho ;
rho_int_trapz = cumtrapz(omega,j_omega_eff ) ;
rho_int = @(x)spline(omega,rho_int_trapz,x).*(x>0) ;
rho = @(x) spline(omega,j_omega_eff,x) ;
rho = @(omega) norm_rho*(sum((0.5.*lambda_D) .* omega_D ./ (omega_D.^2 + omega.^2) ,1) + sum((0.5.*lambda_UBO) .* Omega_UBO.^2 .* gamma_UBO ./ ( (omega.^2-Omega_UBO.^2).^2 + gamma_UBO.^2 .* omega.^2),1));

beta = 1;
N_disc = 1000 ;
omega_max = 50000 ;
N_omega = 160000 ;

[spec,omega_spec,omega_disc,weights,t,c_t_0,g_t] = calculateSpectrumMidPoint(rho_int,Lambda,0,N_disc,N_omega,omega_max,beta,Omega_max) ;
% n_t = N_omega ;
% t_max = 2.0*pi/(2*omega_max) ;
% t = (t_max)*(0:(n_t-1))' ;
% omegas = omega ;
% omega_disc = omegas ;
% weights = j_omega_eff ;
% g_t = -sum(weights ./ omegas ...
%         .*( coth( 0.5*beta * omegas).*(1.-cos(omegas.*t)) + 1.0i * sin(omegas.*t)),2) ;
% c_t_0 = exp( Lambda*g_t ) ;

kappas = [] ;
epsilons = 0:5:300 ;
for i = 1:numel(epsilons)
kappas = [kappas,2*real(trapz(t,c_t_0.*exp(1.0i*epsilons(i)*t))) ] ;
end
% omega_disc = omega(1:end) ; weights = j_omega_eff(1:end) ;

kappa_FGR = 2*real(trapz(t,c_t_0.*exp(1.0i*epsilon*t))) ;
kappa_MT = 2*(sqrt(pi*beta)/sqrt(4*Lambda)) * exp(-(epsilon-Lambda)^2 / (4*Lambda/beta)) ;

k_max = 200 ;
M_max = 1 ;

[nus_D,cs_D,cbars_D] = generateNusAndCsDebye(omega_D,lambda_D,beta,k_max) ;
[nus_UBO,cs_UBO,cbars_UBO] = generateNusAndCsUBO(gamma_UBO,lambda_UBO,Omega_UBO,beta,k_max) ;
% [nus_UBO,cs_UBO] = generateNusAndCsOBO(gamma_OBO,lambda_OBO,Omega_OBO,beta,k_max) ;
% cbars_UBO = conj(cs_UBO) ;
% cs_D = zeros([0,k_max+1]) ; cbars_D= zeros([0,k_max+1]) ; nus_D = zeros([0,k_max+1]) ;
n_UBO = numel(lambda_UBO) ;
n_D = numel(lambda_D) ;
nus = [transpose(nus_D(:,1)),(reshape(transpose(nus_UBO(:,1:2)),[1,2*n_UBO]))] ;
cs = [transpose(cs_D(:,1)),(reshape(transpose(cs_UBO(:,1:2)),[1,2*n_UBO]))] ; 
cbars = [transpose(cs_D(:,1)),(reshape(transpose(cbars_UBO(:,1:2)),[1,2*n_UBO]))] ; 

nus = [nus,nus_UBO(1,(3:(M_max+2)))] ;
cs = [cs ,sum( cs_D(:,2:(M_max+1)),1)+sum( cs_UBO(:,3:(M_max+2),1) )] ;
cbars = [cbars ,sum(cbars_D(:,2:(M_max+1),1))+sum(cbars_UBO(:,3:(M_max+2)),1)] ;


nus_term = nus_UBO(1,(M_max+3):(2+k_max)) ;
cs_term = sum(cs_UBO(:,(M_max+3):(2+k_max)),1) + sum(cs_D(:,(M_max+2):(1+k_max)),1);
cbars_term = sum(cbars_UBO(:,(M_max+3):(2+k_max)),1) + sum(cbars_D(:,(M_max+2):(1+k_max)),1);

Lambda_var = Lambda ;
% epsilon = 20 ;

q = (epsilon-Lambda_var)/(sqrt(2*Lambda_var/beta)) ;
q2 = (-epsilon-Lambda_var)/(sqrt(2*Lambda_var/beta)) ;
Delta_class = sqrt(2*Lambda_var/beta) ;
Delta_QM = sqrt(2*0.5*Lambda_var*sum(weights.*omega_disc.*(coth(beta*omega_disc/2)))) ;
% Delta = Delta_QM ;
% Delta_K = sqrt(real(sum(c_ks .* (( 1-exp( +1.0i*beta*nu_ks ) )./( -1.0i*beta*nu_ks )) ))) ;
% Delta_QM = Delta ;

% Delta_QM = sqrt(2*0.5*Lambda_var*sum(weights.*(cosh(beta*omega_disc/2).^2)).*omega_disc) /2 ;
% Delta_QM = Delta ;
q = (epsilon-Lambda_var)/(Delta_QM) ;
q2 = (-epsilon-Lambda_var)/(Delta_QM) ;
% M = mean(C_ts,2);
% M = M / M(1) ;
% M_class = M ;
% c_0_ps = 299792458e2*1e-12 ; % cm/ps
t_conv = 1 ;
t_old = t ;
t = t' ;

t = linspace(0,50,20000) ;
% C_K = sum(weights.*(cosh(omega_disc*beta)/beta).*cos(omega_disc.*t'*t_conv)./omega_disc,2) ;
% C_K = sum(weights.*(cosh(omega_disc*beta)/beta).*omega_disc.*cos(omega_disc.*t'*t_conv),2) ;
% C_K = sum(weights.*(coth(omega_disc*beta/2)*beta).*cos(omega_disc.*t'*t_conv).*omega_disc,2) ; % NEED TO CHECK THIS real part of CORR not Kubo
% C_K = sum(weights.*(coth(omega_disc*beta/2)*beta).*cos(omega_disc.*t'*t_conv).*omega_disc,2) ; % NEED TO CHECK THIS
% C_K = sum(weights.*(coth(omega_disc*beta/2)*beta).*cos(omega_disc.*t'*t_conv).*sinh(omega_disc*beta),2) ;
% C_K = sum(weights.*(cosh(omega_disc*beta/2).^2*beta).*cos(omega_disc.*t'*t_conv),2) ;
% C_K = sum(weights.*(1*cos(omega_disc.*t'*t_conv)),2) ; % classical limit beta -> 0
% C_K_old = C_K ;
c_ks = [cs,cs_term] ;
nu_ks = [nus,nus_term] ;
% c_ks = cs ;
% nu_ks = nus ;
C_K = sum(c_ks .* (( 1-exp( +1.0i*beta*nu_ks ) )./( -1.0i*beta*nu_ks )) .*exp(-nu_ks.*t'),2) ;
% C_K = sum(c_ks .* (1) .*exp(-nu_ks.*t'),2) ;
% C_K = sum([cs,cs_term].*exp(-[nus,nus_term].*t'),2) ;
C_K = real(C_K) ;



M = (C_K/C_K(1))*(Delta_class/Delta_QM) ;
% M = M' .*(exp(-(t/10).^2)) ;
M = M' .*(0.5 + 0.5*cos((pi)*t/t(end))) ;
% gamma = -(M(2)-M(1))/(dt) ;

f = (1./(sqrt(1-M.^2))).*exp(q^2 * M ./(1+M)) - 1 ;
f2 = (1./(sqrt(1-M.^2))).*exp(q2^2 * M ./(1+M)) - 1 ;

% tau = trapz(t(2:end),f(2:end)) ;
% tau  = tau + exp(0.5*q^2) * (1/sqrt(2*gamma)) * 2 * sqrt(dt) - dt ;
tau = trapz(t(1:end),f(1:end)) ;
% tau  = tau + exp(0.5*q^2) * (1/sqrt(2*gamma)) * 2 * sqrt(2*dt) - 2*dt ;
tau = tau * exp(-0.5*q^2) ;

tau2 = trapz(t(1:end),f2(1:end)) ;
% tau2  = tau2 + exp(0.5*q2^2) * (1/sqrt(2*gamma)) * 2 * sqrt(2*dt) - 2*dt ;
tau2 = tau2 * exp(-0.5*q2^2) ;

tau_tot = tau + tau2 ;
corr = tau_tot * sqrt(2*pi)/Delta_QM ;
corr_0 = corr ;
Delta_QM_0 = Delta_QM ;

k_1 = (2*pi/sqrt(2*pi*Delta_QM^2)) *exp(-0.5*q^2) ;
k_2 = (2*pi/sqrt(2*pi*Delta_QM^2)) *exp(-0.5*q2^2) ;

corr2 = (tau * real(kappa_FGR)/k_1  + tau2 * real(kappa_FGR) * exp(-beta*epsilon)/k_2) * sqrt(2*pi)/Delta_QM ;
% c_0_ps = 299792458e2*1e-12 ; % cm/ps
% hbar = 6.582119569e-16 ; 

% corr = sqrt(2*pi)*(((V)^2)/(Delta*eV_in_cminv*hbar))*tau_s ;
% corr2 = sqrt(2*pi)*(((V)^2)/(Delta*eV_in_cminv*hbar))*tau2_s ;
% 
% % tau_int = trapz(t,M) ;
% % tau_int = (1/13.5175) ;
% tau_int = 0.0660 ;
% M_exact = M ;
% M = exp(-t'/tau_int) ;
% % M = 0.7 * (Lambda_env/Lambda_var) * exp(-t'/tau_int) ;
% f_approx = (1./(sqrt(1-M.^2))).*exp(q^2 * M ./(1+M)) - 1 ;
% tau_approx = exp(-0.5*q^2) * trapz(t,f_approx) ;
% corr_approx = 2*pi*V^2/(abs(epsilon-Lambda_var)*eV_in_cminv*hbar) * tau_int*1e-9 ;
% 
% 
% Delta_QM = sqrt(2.0*0.5*Lambda_var*sum(weights.*omega_disc.*(coth(beta*omega_disc/2)))) ;
% 
k_pade = real(Deltas.^2 * real(kappa_FGR) ./ (1+Deltas.^2 *corr)) ;
k_pade2 = real(Deltas.^2 * real(kappa_FGR) ./ (1+Deltas.^2 *corr2)) ;
k_FGR = Deltas.^2 * real(kappa_FGR) ;

k_epscorr = [] ;
k_epscorr2 = [] ;
k_epscorr3 = [] ;
k_epscorrpade = [] ;
k_epscorrpade2 = [] ;
j_eff_A = imag(trapz(t_old,(c_t_0).*exp(1.0i*epsilon*t_old))) ;
j_eff_B = imag(trapz(t_old,conj(c_t_0).*exp(1.0i*epsilon*t_old))) ;
eta_mins = [] ;
% eta_min = -0.5 ;
eta_min = -1 + Deltas(end)/(2*Lambda+0*Deltas(end)) ;
sin_2theta = sin(2*0.5*atan(2*Deltas(end)/(epsilon-0.5*Lambda))) ;
theta_opts = [] ;
theta_opt = 1e-3 ;
% for n = numel(Deltas):-1:1
    for n = 1:numel(Deltas)
    % epsilon_eff = epsilon + Deltas(n)^2/(epsilon-Lambda) + Deltas(n)^2/(epsilon+Lambda) ;
    epsilon_eff = epsilon + Deltas(n)^2 * (j_eff_A + j_eff_B) ;
    epsilon_eff = epsilon + sqrt(Deltas(n)^2 + 0.25*(epsilon-Lambda)^2)-0.5*(epsilon-Lambda) + sqrt(Deltas(n)^2 + 0.25*(epsilon+Lambda)^2)-0.5*(epsilon+Lambda) ;
    Lambda_eff = Lambda - (sqrt(Deltas(n)^2 + 0.25*(epsilon-Lambda)^2)-0.5*(epsilon-Lambda)) + sqrt(Deltas(n)^2 + 0.25*(epsilon+Lambda)^2)-0.5*(epsilon+Lambda) ;
    
    % epsilon_eff = epsilon + Deltas(n)^2 * (j_eff_A + j_eff_B) ;
    % Lambda_eff = Lambda + Deltas(n)^2 * (-j_eff_A + j_eff_B)  ;
    
    epsilon_eff = sqrt(epsilon^2 + 4 * Deltas(n)^2) ;
    Lambda_eff =  Lambda / sqrt(1 + 4 * (Deltas(n)/epsilon)^2) ;
    
    Delta = Deltas(n) ;
    Lambda_theta = @(theta) (cos(2*theta).^2) * Lambda ;
    epsilon_theta = @(theta)  cos(2*theta)*epsilon + 2*sin(2*theta)*Delta ;
    % k_theta0 = @(theta) (sec(2*theta).^2).*(2*Delta^2).*real(trapz(t_old,exp( Lambda_theta(theta).*g_t ).*exp(1.0i*epsilon_theta(theta).*t_old).*(0.5 + 0.5*cos((pi)*t_old/t_old(end)))) );
    % k_theta = @(theta) (k_theta0(theta)>0).*k_theta0(theta)>0 + (k_theta0(theta)<0) * 1e4 ;
    k_theta_f = @(theta) (sec(2*theta).^2).*(2*Delta^2).*real(trapz(t_old,exp( Lambda_theta(theta).*g_t ).*exp(1.0i*epsilon_theta(theta).*t_old).*(0.5 + 0.5*cos(0*(pi)*t_old/t_old(end)))) );
    % k_theta = @(theta) (sec(2*theta).^2).*(2*Delta^2).*real(trapz(t_old,exp( Lambda_theta(theta).*g_t ).*exp(1.0i*epsilon_theta(theta).*t_old))) ;
    
    theta_eta = @(eta)0.5 * atan(2*Delta./(epsilon+eta*Lambda)) ;
    
    % eta_0 = eta_min ;
    % eta_0 = 0 ;
    % opts = optimoptions("fmincon","Display","off") ;
    % opts = optimoptions("patternsearch","FunctionTolerance",1e-8,"MeshTolerance",1e-4) ;
    % % [eta_min,k_theta_min] = patternsearch(@(eta)k_theta(theta_eta(eta)),eta_0,[],[],[],[],[-0.99],[2],[],opts)  ;
    % theta_0 = 0.5 * atan(2*Deltas(n)/(epsilon-0.99*Lambda)) ;
    % % theta_0 = pi/4-1e-2 ;
    % [sin_2theta,k_theta_min] = patternsearch(@(s)(k_theta(asin(s)/2)),sin_2theta,[],[],[],[],[0],[sin(2*theta_0)],[],opts)  ;
    % [sin_2theta0,k_theta_min0] = patternsearch(@(s)(k_theta(asin(s)/2)),0,[],[],[],[],[0],[sin(2*theta_0)],[],opts)  ;
    % if k_theta_min0 < k_theta_min 
    %     sin_2theta = sin_2theta0 ;
    % end
    % theta = asin(sin_2theta)/2 ;
    
    % etas = [linspace(-1.0,-0.5,20),linspace(-0.5,2,80)] ;
    % k_thetas = k_theta(theta_eta(etas)) ;
    % inds = k_thetas > 0;
    % k_thetas = k_thetas(inds) ; etas = etas(inds) ;
    % [k_min,i] = min(k_thetas) ;
    % f = @(eta) spline(etas,k_thetas,eta) ;
    % eta_min = etas(i) ;
    % % [eta_min] = fmincon(f,etas(i),[],[],[],[],[-1],[2],[],opts)  ;
    % theta = theta_eta(eta_min) ;
    % eta_mins(n) = eta_min ;

    % thetas = linspace(0,0.5*atan(2*Delta/(epsilon-Lambda+1e-4)),30) ;
    thetas = linspace(0.5*atan(2*Delta/(epsilon+Lambda)),0.5*atan(2*Delta/(epsilon-Lambda+1e-4)),8) ;
   
    % thetas = linspace(-pi/4,pi/4,60) ;
    k_thetas = k_theta_f(thetas) ;
    if abs(Lambda-epsilon) < 0.01 
        k_thetas(end) = (2*pi)*0.25*Delta*Lambda*rho(2*Delta)*coth(beta*(2*Delta)/2) ;
    end
    inds = k_thetas > 0;
    k_thetas = k_thetas(inds) ; thetas = thetas(inds) ;
    
    [k_min,i] = min(k_thetas) ;
    % f = @(eta) spline(etas,k_thetas,eta) ;
    theta = thetas(i) ;
    
    
    % V_A = @(x)(1/(4*Lambda))*(x-(epsilon+Lambda)).^2 ;
    % V_B = @(x)(1/(4*Lambda))*(x-(epsilon-Lambda)).^2 +epsilon ;
    % V_plus = @(x) (V_A(x)+V_B(x))/2 + 0.5*sqrt((2*Delta)^2 + x.^2) ;
    % opts = optimoptions("fminunc","OptimalityTolerance",1e-16) ;
    % x_min = fminunc(V_plus,0.1,opts) ;
    % theta = 0.5 * atan(2*Deltas(n)/(x_min)) ;
    % theta = 0.5 * atan(2*Deltas(n)/(epsilon-1*Lambda)) ;
    % theta = 0.5 * atan(2*Deltas(n)/(epsilon-0.5*Lambda)) ;
    % theta = 0.5 * atan(2*Deltas(n)*j_eff_A) ;
    % if abs(epsilon-Lambda) < 1e-10 
    %     theta = 0.5 * atan(inf) ;
    % end
    % theta = 0.5*atan((2*Delta)^10 / ((epsilon-Lambda)^10 + (epsilon*Delta)^5)) ;
    
    % f = @(x) 2* Delta * sec(2*x)./(epsilon * cos(2*x) + Delta * sin(2*x) - Lambda * (cos(2*x)).^2) +0.5*tan(2*x)*(epsilon * cos(2*x) + Delta * sin(2*x) - Lambda * (cos(2*x)).^2) - tan(2*x) ;
    % theta_opt = fsolve(f,0.5*atan(2*Delta/epsilon)) ;
    % theta = theta_opt ;
    
    % energy gap at + at min
    % f = @(x) 4* Delta * cos(2*x) - 2*epsilon * sin(2*x) + 4* Lambda * sin(2*x).* cos(2*x) ;
    % absolute energy of + at min
    % f = @(x) -epsilon * sin(2*x) + cos(2*x).*(2*Delta + Lambda * sin(2*x)) ;
     % f = @(x) -epsilon * sin(2*x) + cos(2*x).*(2*Delta + Lambda * sin(2*x)) ;
    % f = @(x) -4* Delta * cos(2*x) + 2* epsilon * sin(2*x) +  Lambda * sin(4*x);
    % f = @(x) 2* Delta * sin(2*x) +  epsilon * cos(2*x) +  Lambda * cos(4*x);
    % f = @(x) 4* Delta * cos(2*x) -2* epsilon * sin(2*x) +  2*Lambda * sin(4*x);
    % f = @(x) -2* Delta * cos(2*x) + epsilon * sin(2*x) -  0.5*Lambda * sin(4*x);
    % theta = fsolve(f,0.5*atan(2*Delta/(epsilon-Lambda)) ) ;
    % if theta >(pi/4) 
    %     theta = fsolve(f,0) ;
    % end

    % theta = 0.5*atan(2*Delta/(epsilon-Lambda)) ;
    % theta = 0.5*atan(2*Delta*j_eff_A) ;
    % theta = 0.5*atan(2*Delta/epsilon) ;
    epsilon_eff = cos(2*theta)*epsilon + 2*sin(2*theta)*Deltas(n); 
    Lambda_eff = (cos(2*theta)^2) * Lambda ;
    % epsilon_theta = @(theta)  cos(2*theta)*epsilon + 2*sin(2*theta)*Delta ;
    
    
    k_epscorr(n) = Deltas(n).^2 *2*real(trapz(t_old,c_t_0.*exp(1.0i*epsilon_eff*t_old) .*(0.5 + 0.5*cos(0*(pi)*t_old/t_old(end))) )) ;
    
    c_t_eff = exp( Lambda_eff*g_t ) ;
    k_epscorr(n) = Deltas(n).^2 *2*real(trapz(t_old,c_t_eff.*exp(1.0i*epsilon_eff*t_old) .*(0.5 + 0.5*cos(0*(pi)*t_old/t_old(end))) )) ;
    k_epscorr2(n) = k_epscorr(n) * (sec(2*theta)^2) ;
    k_epscorr2(n) = k_min ;
    DA_eff = epsilon + sqrt(Deltas(n)^2 + 0.25*(epsilon-Lambda)^2)-0.5*(epsilon-Lambda) + sqrt(Deltas(n)^2 + 0.25*(epsilon+Lambda)^2)-0.5*(epsilon+Lambda) ;
    
    k_epscorr3(n) = k_epscorr2(n) * (1+exp(-beta*epsilon_eff)) / (1+exp(-beta*DA_eff)) ;
    
    Lambda_var = Lambda_eff ;
    Delta_QM = sqrt(Lambda_eff/Lambda) * Delta_QM_0 ;
    
    q = (epsilon_eff-Lambda_var)/(Delta_QM) ;
    q2 = (-epsilon_eff-Lambda_var)/(Delta_QM) ;
    
    f = (1./(sqrt(1-M.^2))).*exp(q^2 * M ./(1+M)-0.5*q^2) - exp(-0.5*q^2) ;
    f2 = (1./(sqrt(1-M.^2))).*exp(q2^2 * M ./(1+M)-0.5*q2^2) - exp(-0.5*q2^2) ;
    % hold on
    % plot(t,f) 
    
    tau = trapz(t(1:end),f(1:end)) ;
    % tau = tau * exp(-0.5*q^2) ;
    
    tau2 = trapz(t(1:end),f2(1:end)) ;
    % tau2 = tau2 * exp(-0.5*q2^2) ;
    
    tau_tot = tau + tau2 ;
    corr = tau_tot * sqrt(2*pi)/Delta_QM  ;
    % corr2 = corr  * (sqrt(1 + 4 * (Deltas(n)/epsilon)^2) + 0*(epsilon_eff)/epsilon)^2 ;
    corr2 = corr * ((sec(2*theta) )^2) ;
    
    k_epscorrpade(n) = k_epscorr(n)* (sec(2*theta) )^2 / (1 + (Deltas(n)^2) * corr) ;
    % gamma = 1 - k_epscorr(n) / ((sqrt(2*pi) * Deltas(n)^2 * exp(-q^2/2))+k_epscorr(n)) ;
    % k_epscorrpade2(n) = k_epscorr(n) * ( 1-gamma + ((gamma) / (1 + (Deltas(n)^2) * corr)) )  ;
    k_epscorrpade2(n) = k_epscorr(n)* (sec(2*theta) )^2 / (1 + (Deltas(n)^2) * corr2)  ;
    theta_opts(n) = theta ;
end
