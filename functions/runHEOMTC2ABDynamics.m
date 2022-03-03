function [O_t,t,junk] = runHEOMTC2ABDynamics(full_system,heom_dynamics,AB_coupling_info)
% convert the bath info into a more use-able form - these are the
% explicitly treated baths
[heom_bath_info_A,heom_bath_info_B] = getBathInformationAB(full_system) ;
% construct the HEOM dynamics genrators for the A & B spaces
[L_A,heom_structure_A] = constructHEOMGenerator(full_system.H_sys_A,heom_bath_info_A, ...
    heom_dynamics.heom_truncation) ;
[L_B,heom_structure_B] = constructHEOMGenerator(full_system.H_sys_B,heom_bath_info_B, ...
    heom_dynamics.heom_truncation, heom_structure_A) ;

% get the dimensions of things
d_heom_A = size(L_A,1) ;
d_heom_B = size(L_B,1) ;
d_hilb_A = size(full_system.H_sys_A,1) ;
d_liou_A = d_hilb_A * d_hilb_A ;
d_hilb_B = size(full_system.H_sys_B,1) ;
d_liou_B = d_hilb_B * d_hilb_B ;
d_heom = d_heom_A + d_heom_B ;

L = sparse([],[],[],d_heom,d_heom) ;
A_inds = 1:d_heom_A ;
B_inds = (d_heom_A+1):d_heom ;
L(A_inds,A_inds) = L_A ;
L(B_inds,B_inds) = L_B ;

[K,c_ts,ts] = constructTC2ABOperator(AB_coupling_info, full_system.H_sys_A, full_system.H_sys_B,d_heom_A,d_heom_B,...
    heom_structure_A,heom_structure_B,heom_bath_info_A.Vs,heom_bath_info_B.Vs,heom_bath_info_A,heom_dynamics.heom_truncation) ;
if (AB_coupling_info.method == "full NZ")
    L = full(L+K) ;
else
    L = L + K ;
end

% construct the rho_0 for the full hierarchy
rho_0_heom = zeros([d_heom,1]) ;
rho_0_heom(1:d_liou_A) = convertToLiouvilleVector(heom_dynamics.rho_0_sys_A) ;
rho_0_heom((d_heom_A+1):(d_heom_A+d_liou_B)) = convertToLiouvilleVector(heom_dynamics.rho_0_sys_B) ;

% set up the observable operators
n_obs_A = numel(heom_dynamics.observables.system_A) ;
n_obs_B = numel(heom_dynamics.observables.system_B) ;
n_obs = n_obs_A + n_obs_B ;
O = sparse([],[],[],n_obs,d_heom) ;
for n = 1:n_obs_A
    O(n,1:d_liou_A) = convertToLiouvilleVector(heom_dynamics.observables.system_A{n})' ;
end
for n = 1:n_obs_B
    O(n+n_obs_A,(d_heom_A+1):(d_heom_A+d_liou_B)) = convertToLiouvilleVector(heom_dynamics.observables.system_B{n})' ;
end
% run the dynamics
integrator = heom_dynamics.integrator ;
if (heom_dynamics.integrator.method == 'SIA')
    [O_t,t] = runDynamicsSIADensityOperator(rho_0_heom,L,...
        integrator.n_steps,integrator.dt,O,integrator.krylov_dim,...
        integrator.krylov_tol) ;
end
junk = K ;
end