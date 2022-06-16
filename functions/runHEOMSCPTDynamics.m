function [O_t,t,L,junk] = runHEOMSCPTDynamics(full_system,heom_dynamics)
% convert the bath info into a more use-able form - these are the
% explicitly treated baths
heom_bath_info_blocks = getBathInformationSCPT(full_system) ;
% get the number of SC blocks
n_blocks = numel(full_system.H_sys) ;
% construct the HEOM dynamics genrators for the different blocks
L_blocks = cell([n_blocks,1]) ; heom_structure_blocks = cell([n_blocks,1]) ;
for k = 1:1
    [L_block,heom_structure_block] = constructHEOMGenerator(full_system.H_sys{k},heom_bath_info_blocks{k}, ...
        heom_dynamics.heom_truncation) ;
    L_blocks{k} = L_block ;
    heom_structure_blocks{k} = heom_structure_block ;
end
for k = 2:n_blocks
    [L_block,heom_structure_block] = constructHEOMGenerator(full_system.H_sys{k},heom_bath_info_blocks{k}, ...
        heom_dynamics.heom_truncation,heom_structure_blocks{1}) ;
    L_blocks{k} = L_block ;
    heom_structure_blocks{k} = heom_structure_block ;
end
% construct the diagonal elements of the generator for the HEOMSCPT object
n_ados = [] ;
d_heoms = [] ;
d_lious = [] ;
block_indices = {};
for k = 1:n_blocks
    n_ados = [n_ados,heom_structure_blocks{k}.n_ados] ;
    d_heoms = [d_heoms,heom_structure_blocks{k}.d_heom] ;
    d_lious = [d_lious,heom_structure_blocks{k}.d_liou] ;
    block_indices{k} = (1+sum(d_heoms(1:end-1))):sum(d_heoms(1:end)) ;
end
d = sum(d_heoms) ;
L = sparse([],[],[],d,d) ;
for k = 1:n_blocks
    L(block_indices{k},block_indices{k}) = L_blocks{k} ;
end

% next add the markovian rate matrix
[K,scpt_junk] = constructHEOMSCPTCouplingOperator(full_system,heom_structure_blocks,heom_bath_info_blocks,heom_dynamics.block_coupling,heom_dynamics.heom_truncation) ;
L = L + K ;
drawnow ;

% add the incoherent processes term if it is specfiied
if isfield(full_system,'incoh_processes')
    K_incoh = constructIncoherentRateOperator(full_system.incoh_processes,n_ados,d_lious) ;
    L = L + K_incoh ;
end


% construct the rho_0 for the full hierarchy
rho_0_heom = zeros([d,1]) ;
d_start = 0 ;
for k = 1:n_blocks
    inds = (d_start+1):(d_start+d_lious(k)) ;
    rho_0_heom(inds) = convertToLiouvilleVector(heom_dynamics.rho_0_sys{k}) ;
    d_start = d_start + d_heoms(k) ;
end

% set up the observable operators
n_obs_blocks = zeros([n_blocks,1]) ;
for k = 1:n_blocks
    n_obs_blocks(k) = numel(heom_dynamics.observables.block{k});
end
n_obs = sum(n_obs_blocks) ;
O = sparse([],[],[],n_obs,d) ;
for k = 1:n_blocks
    n_start = sum(n_obs_blocks(1:(k-1))) ;
    d_start = sum(d_heoms(1:(k-1))) ;
    for n = 1:n_obs_blocks(k)
        O(n+n_start,(1+d_start):(d_lious(k)+d_start)) = convertToLiouvilleVector(heom_dynamics.observables.block{k}{n})' ;
    end
end

% run the dynamics
integrator = heom_dynamics.integrator ;
fprintf("Starting dynamics.\n");
if (heom_dynamics.integrator.method == 'SIA')
    [O_t,t, rho_final_heom] = runDynamicsSIADensityOperator(rho_0_heom,L,...
        integrator.n_steps,integrator.dt,O,integrator.krylov_dim,...
        integrator.krylov_tol) ;
elseif (integrator.method == 'adaptive taylor')
    [O_t, t, rho_final_heom] = runDynamicsAdaptiveTaylorDensityOperator(rho_0_heom,L,...
        integrator.t_max,O,integrator.order,integrator.tol) ;
elseif (integrator.method == 'adaptive SIA')
    [O_t,t,rho_final_heom] = runDynamicsSIAAdaptiveDensityOperator(rho_0_heom,L,...
        integrator.n_steps,integrator.dt,O,integrator.krylov_dim,...
        integrator.krylov_tol) ;
end

junk = scpt_junk ;
end