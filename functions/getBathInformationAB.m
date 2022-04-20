function [heom_bath_info_A,heom_bath_info_B,lambda_Ds,omega_Ds,lambda_OBOs,Omega_OBOs,gamma_OBOs,lambda_UBOs,Omega_UBOs,gamma_UBOs,V_As,V_Bs] = getBathInformation(full_system)
% this method reorganises information on the baths into a form used to
% actually construct the HEOM, returned in the bath_info struct
baths = full_system.baths ;
n_baths = numel(baths) ;

lambda_Ds = [] ; omega_Ds = [] ;
lambda_OBOs = [] ; Omega_OBOs = [] ; gamma_OBOs = [] ;
lambda_UBOs = [] ; Omega_UBOs = [] ; gamma_UBOs = [] ;
V_As = {} ;
V_Bs = {} ;


for i = 1:n_baths
    if (baths{i}.spectral_density == "debye")
        lambda_Ds = [lambda_Ds,baths{i}.lambda_D] ;
        omega_Ds = [omega_Ds,baths{i}.omega_D] ;
        V_As = [V_As,{baths{i}.V_A}] ; 
        V_Bs = [V_Bs,{baths{i}.V_B}] ; 
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "OBO")
        lambda_OBOs = [lambda_OBOs,baths{i}.lambda] ;
        Omega_OBOs = [Omega_OBOs,baths{i}.Omega] ;
        gamma_OBOs = [gamma_OBOs,baths{i}.gamma] ;
        V_As = [V_As,{baths{i}.V_A}] ; 
        V_Bs = [V_Bs,{baths{i}.V_B}] ; 
    end
end

for i = 1:n_baths
    if (baths{i}.spectral_density == "UBO")
        lambda_UBOs = [lambda_UBOs,baths{i}.lambda] ;
        Omega_UBOs = [Omega_UBOs,baths{i}.Omega] ;
        gamma_UBOs = [gamma_UBOs,baths{i}.gamma] ;
        V_As = [V_As,{baths{i}.V_A}] ; 
        V_Bs = [V_Bs,{baths{i}.V_B}] ;  
    end
end

heom_bath_info = struct ;
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

heom_bath_info_A = heom_bath_info ;
heom_bath_info_A.Vs = V_As ;
heom_bath_info_B = heom_bath_info ;
heom_bath_info_B.Vs = V_Bs ;


end