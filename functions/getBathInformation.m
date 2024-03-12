function [heom_bath_info,lambda_Ds,omega_Ds,lambda_OBOs,Omega_OBOs,gamma_OBOs,lambda_UBOs,Omega_UBOs,gamma_UBOs,Vs] = getBathInformation(full_system)
% this method reorganises information on the baths into a form used to
% actually construct the HEOM, returned in the bath_info struct
baths = full_system.baths ;
n_baths = numel(baths) ;

lambda_Ds = [] ; omega_Ds = [] ;
lambda_OBOs = [] ; Omega_OBOs = [] ; gamma_OBOs = [] ;
lambda_UBOs = [] ; Omega_UBOs = [] ; gamma_UBOs = [] ;
lambda_Ds_pade = [] ; omega_Ds_pade = [] ; pade_approximants = [] ; N_pade = [] ;
Vs = {} ;
nus_custom = {} ; cs_custom = {} ; cbars_custom = {} ; cs_trunc_custom = {} ; nus_trunc_custom = {} ; cbars_trunc_custom = {} ;

for i = 1:n_baths
    if (baths{i}.spectral_density == "debye")
        lambda_Ds = [lambda_Ds,baths{i}.lambda_D] ;
        omega_Ds = [omega_Ds,baths{i}.omega_D] ;
        Vs = [Vs,{baths{i}.V}] ; 
    end
end


for i = 1:n_baths
    if (baths{i}.spectral_density == "OBO")
        lambda_OBOs = [lambda_OBOs,baths{i}.lambda] ;
        Omega_OBOs = [Omega_OBOs,baths{i}.Omega] ;
        gamma_OBOs = [gamma_OBOs,baths{i}.gamma] ;
        Vs = [Vs,{baths{i}.V}] ; 
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "UBO")
        lambda_UBOs = [lambda_UBOs,baths{i}.lambda] ;
        Omega_UBOs = [Omega_UBOs,baths{i}.Omega] ;
        gamma_UBOs = [gamma_UBOs,baths{i}.gamma] ;
        Vs = [Vs,{baths{i}.V}] ; 
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

for i = 1:n_baths
    if (baths{i}.spectral_density == "custom")
        nus_custom = [nus_custom,{baths{i}.nus}] ; 
        cs_custom = [cs_custom,{baths{i}.cs}] ; 
        cbars_custom = [cbars_custom,{baths{i}.cbars}] ; 
        cs_trunc_custom = [cs_trunc_custom,{baths{i}.cs_trunc}] ; 
        cbars_trunc_custom = [cbars_trunc_custom,{baths{i}.cbars_trunc}] ; 
        nus_trunc_custom = [nus_trunc_custom,{baths{i}.nus_trunc}] ; 
        Vs = [Vs,{baths{i}.V}] ; 
    end
end

heom_bath_info = struct() ;
heom_bath_info.n_baths = n_baths ;
heom_bath_info.Vs = Vs ;
heom_bath_info.lambda_Ds = lambda_Ds ;
heom_bath_info.omega_Ds = omega_Ds ;
heom_bath_info.lambda_OBOs = lambda_OBOs ;
heom_bath_info.Omega_OBOs = Omega_OBOs ;
heom_bath_info.gamma_OBOs = gamma_OBOs ;
heom_bath_info.lambda_UBOs = lambda_UBOs ;
heom_bath_info.Omega_UBOs = Omega_UBOs ;
heom_bath_info.gamma_UBOs = gamma_UBOs ;
heom_bath_info.beta = full_system.beta ;
heom_bath_info.lambda_Ds_pade = lambda_Ds_pade ;
heom_bath_info.omega_Ds_pade = omega_Ds_pade ;
heom_bath_info.N_pade = N_pade ;
heom_bath_info.pade_approximants = pade_approximants ;
heom_bath_info.nus_custom = nus_custom ;
heom_bath_info.cs_custom = cs_custom ;
heom_bath_info.cbars_custom = cbars_custom ;
heom_bath_info.nus_trunc_custom = nus_trunc_custom ;
heom_bath_info.cs_trunc_custom = cs_trunc_custom ;
heom_bath_info.cbars_trunc_custom = cbars_trunc_custom ;


end