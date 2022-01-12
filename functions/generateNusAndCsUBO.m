function [nus,cs,cbars] = generateNusAndCsUBO(gamma_Bs,lambda_Bs,Omega_Bs,beta,M)
% generates the nu_jks and c_jks for a set of Debye baths with debye
% frequencies specified by the row vector gamma_Ds and temperature beta

% get the number of baths
N_baths = length(gamma_Bs) ;
% empty arrays for the c_jk and nu_jk values
nus = zeros(N_baths,(M+2)) ;
cs = zeros(N_baths,(M+2)) ; 
cbars = zeros(N_baths,(M+2)) ; 

lambda_Bs_col = reshape(lambda_Bs,[N_baths,1]) ;
gamma_Bs_col = reshape(gamma_Bs,[N_baths,1]) ; 
Omega_Bs_col = reshape(Omega_Bs,[N_baths,1]) ; 
kappas_col = sqrt(-(0.5 * gamma_Bs_col).^2 + Omega_Bs_col.^2) ;

% create the nu_jks according to https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.104.250401
nus(:,1) = 0.5*gamma_Bs_col + 1.0i * kappas_col ; % nu_+
nus(:,2) = 0.5*gamma_Bs_col - 1.0i * kappas_col ; % nu_-
nus(:,3:(M+2)) = (2*pi/beta) * repmat((1:M),[N_baths,1]) ;

% create the cs according to the above()

cs(:,1) = -(0.5*lambda_Bs_col.* (Omega_Bs_col.^2)./kappas_col).*(coth(0.5i*beta*nus(:,1))-1.0) ;
cs(:,2) = (0.5*lambda_Bs_col.* (Omega_Bs_col.^2)./kappas_col).*(coth(0.5i*beta*nus(:,2))-1.0)  ;
cs(:,3:(M+2)) = (4.0/beta)*(gamma_Bs_col.*lambda_Bs_col.*(Omega_Bs_col.^2)).*nus(:,3:(M+2))./((gamma_Bs_col.^2).*(nus(:,3:(M+2)).^2) - (nus(:,3:(M+2)).^2 + Omega_Bs_col.^2).^2) ;

% also need an array of cbars for the UBO
cbars = cs ;
cbars(:,1) = cs(:,2) ;
cbars(:,2) = cs(:,1) ;

end