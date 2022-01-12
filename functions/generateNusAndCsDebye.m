function [nus,cs] = generateNusAndCsDebye(gamma_Ds,lambda_Ds,beta,M)
% generates the nu_jks and c_jks for a set of Debye baths with debye
% frequencies specified by the row vector gamma_Ds and temperature beta

% get the number of baths
N_baths = length(gamma_Ds) ;
% empty arrays for the c_jk and nu_jk values
nus = zeros(N_baths,(M+1)) ;
cs = zeros(N_baths,(M+1)) ; 

% create the nu_jks according to https://journals.aps.org/prl/pdf/10.1103/PhysRevLett.104.250401
nus(:,1) = reshape(gamma_Ds,[N_baths,1]) ;
nus(:,2:(M+1)) = (2*pi/beta) * repmat((1:M),[N_baths,1]) ;

% create the cs according to the above()
lambda_Ds_col = reshape(lambda_Ds,[N_baths,1]) ;
gamma_Ds_col = reshape(gamma_Ds,[N_baths,1]) ; 
cs(:,1) = gamma_Ds_col.*lambda_Ds_col .* (-1.0i + cot(0.5*beta*gamma_Ds_col)) ;
cs(:,2:(M+1)) = (4.0/beta)*gamma_Ds_col.*lambda_Ds_col.*nus(:,2:(M+1))./(nus(:,2:(M+1)).^2 - gamma_Ds_col.^2) ;

end