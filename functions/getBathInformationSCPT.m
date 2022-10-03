function [heom_bath_info_blocks,lambda_Ds,omega_Ds,lambda_OBOs,Omega_OBOs,gamma_OBOs,lambda_UBOs,Omega_UBOs,gamma_UBOs,V_As,V_Bs] = getBathInformationSCPT(full_system)
% this method reorganises information on the baths into a form used to
% actually construct the HEOM, returned in the bath_info struct
baths = full_system.baths ;
n_baths = numel(baths) ;
n_blocks = numel(full_system.H_sys) ;

lambda_Ds = [] ; omega_Ds = [] ;
lambda_OBOs = [] ; Omega_OBOs = [] ; gamma_OBOs = [] ;
lambda_UBOs = [] ; Omega_UBOs = [] ; gamma_UBOs = [] ;
lambda_Ds_pade = [] ; omega_Ds_pade = [] ; pade_approximants = [] ; N_pade = [] ;
Vs = {} ;

for i = 1:n_baths
    if (baths{i}.spectral_density == "debye")
        lambda_Ds = [lambda_Ds,baths{i}.lambda_D] ;
        omega_Ds = [omega_Ds,baths{i}.omega_D] ;
%         Vs = [Vs,{baths{i}.V}] ;
        Vs = [Vs,{full_system.Vs{i}}] ;
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "OBO")
        lambda_OBOs = [lambda_OBOs,baths{i}.lambda] ;
        Omega_OBOs = [Omega_OBOs,baths{i}.Omega] ;
        gamma_OBOs = [gamma_OBOs,baths{i}.gamma] ;
%         Vs = [Vs,{baths{i}.V}] ;
        Vs = [Vs,{full_system.Vs{i}}] ;
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "UBO")
        lambda_UBOs = [lambda_UBOs,baths{i}.lambda] ;
        Omega_UBOs = [Omega_UBOs,baths{i}.Omega] ;
        gamma_UBOs = [gamma_UBOs,baths{i}.gamma] ;
        [Vs,{full_system.Vs{i}}] ;
%         Vs = [Vs,{baths{i}.V}] ;
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "debye (pade)")
        lambda_Ds_pade = [lambda_Ds_pade,baths{i}.lambda_D] ;
        omega_Ds_pade = [omega_Ds_pade,baths{i}.omega_D] ;
        pade_approximants = [pade_approximants,baths{i}.approximant_type ]; 
        N_pade = [N_pade,baths{i}.N_pade] ;
        Vs = [Vs,{baths{i}.V}] ; 
    end
end

heom_bath_info_blocks = {} ;
for k = 1:n_blocks
    heom_bath_info = struct() ;
    heom_bath_info.n_baths = n_baths ;
    heom_bath_info.lambda_Ds = lambda_Ds ;
    heom_bath_info.omega_Ds = omega_Ds ;
    heom_bath_info.lambda_OBOs = lambda_OBOs ;
    heom_bath_info.Omega_OBOs = Omega_OBOs ;
    heom_bath_info.gamma_OBOs = gamma_OBOs ;
    heom_bath_info.lambda_UBOs = lambda_UBOs ;
    heom_bath_info.Omega_UBOs = Omega_UBOs ;
    heom_bath_info.gamma_UBOs = gamma_UBOs ;
    heom_bath_info.beta = full_system.beta ;
    heom_bath_info.Vs = {} ;
    for j = 1:n_baths
        heom_bath_info.Vs{j} = [Vs{j}{k}] ;
    end
    heom_bath_info.lambda_Ds_pade = lambda_Ds_pade ;
    heom_bath_info.omega_Ds_pade = omega_Ds_pade ;
    heom_bath_info.N_pade = N_pade ;
    heom_bath_info.pade_approximants = pade_approximants ;
    heom_bath_info_blocks = [heom_bath_info_blocks,{heom_bath_info}] ;
end


end