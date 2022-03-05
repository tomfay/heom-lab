function [O_t,t,rho_final_heom,L_heom] = runHEOMDynamics(full_system,heom_dynamics)
% convert the bath info into a more use-able form
heom_bath_info = getBathInformation(full_system) ;
% construct the HEOM dynamics genrator as a sparse matrix
[L_heom,ado_indices] = constructHEOMGenerator(full_system.H_sys,heom_bath_info, ...
    heom_dynamics.heom_truncation) ;
% get the dimensions of things
d_heom = size(L_heom,1) ;
d_hilb = size(full_system.H_sys,1) ;
d_liou = d_hilb * d_hilb ;
% construct the rho_0 for the full hierarchy
rho_0_heom = zeros([d_heom,1]) ;
rho_0_heom(1:d_liou) = convertToLiouvilleVector(heom_dynamics.rho_0_sys) ;
% set up the observable operators
n_obs = numel(heom_dynamics.observables.system) ;
O = sparse([],[],[],n_obs,d_heom) ;
for n = 1:n_obs
    O(n,1:d_liou) = convertToLiouvilleVector(heom_dynamics.observables.system{n})' ;
end
% run the dynamics
integrator = heom_dynamics.integrator ;
if (integrator.method == 'SIA')
    [O_t,t,rho_final_heom] = runDynamicsSIADensityOperator(rho_0_heom,L_heom,...
        integrator.n_steps,integrator.dt,O,integrator.krylov_dim,...
        integrator.krylov_tol) ;
elseif (integrator.method == 'adaptive taylor')
    [O_t, t, rho_final_heom] = runDynamicsAdaptiveTaylorDensityOperator(rho_0_heom,L_heom,...
        integrator.t_max,O,integrator.order,integrator.tol) ;
elseif (integrator.method == 'adaptive SIA')
    [O_t,t,rho_final_heom] = runDynamicsSIAAdaptiveDensityOperator(rho_0_heom,L_heom,...
        integrator.n_steps,integrator.dt,O,integrator.krylov_dim,...
        integrator.krylov_tol) ;

end

end