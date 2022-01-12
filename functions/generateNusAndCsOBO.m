function [nus,cs] = generateNusAndCsOBO(gamma_Bs,lambda_Bs,Omega_Bs,beta,M)
% generates the nu_jks and c_jks for a set of Debye baths with debye
% frequencies specified by the row vector gamma_Ds and temperature beta

% get the number of baths
N_baths = length(gamma_Bs) ;
% empty arrays for the c_jk and nu_jk values
nus = zeros(N_baths,(M+2)) ;
cs = zeros(N_baths,(M+2)) ; 

lambda_Bs_col = reshape(lambda_Bs,[N_baths,1]) ;
gamma_Bs_col = reshape(gamma_Bs,[N_baths,1]) ; 
Omega_Bs_col = reshape(Omega_Bs,[N_baths,1]) ; 
kappas_col = sqrt((0.5 * gamma_Bs_col).^2 - Omega_Bs_col.^2) ;

% create the nu_jks according to https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.104.250401
nus(:,1) = 0.5*gamma_Bs_col + kappas_col ;
nus(:,2) = 0.5*gamma_Bs_col - kappas_col ;
nus(:,3:(M+2)) = (2*pi/beta) * repmat((1:M),[N_baths,1]) ;

% create the cs according to the above()

cs(:,1) = -(0.5*lambda_Bs_col.* (Omega_Bs_col.^2)./kappas_col).*(cot(0.5*beta*nus(:,1))-1.0i) ;
cs(:,2) = (0.5*lambda_Bs_col.* (Omega_Bs_col.^2)./kappas_col).*(cot(0.5*beta*nus(:,2))-1.0i)  ;
cs(:,3:(M+2)) = (4.0/beta)*(gamma_Bs_col.*lambda_Bs_col.*(Omega_Bs_col.^2)).*nus(:,3:(M+2))./((gamma_Bs_col.^2).*(nus(:,3:(M+2)).^2) - (nus(:,3:(M+2)).^2 + Omega_Bs_col.^2).^2) ;

end