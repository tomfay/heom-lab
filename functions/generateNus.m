function nus = generateNus(gamma_Ds,beta,M)

N_baths = length(gamma_Ds) ;
nus = zeros(N_baths,(M+1)) ;

nus(:,1) = reshape(gamma_Ds,[N_baths,1]) ;
nus(:,2:(M+1)) = (2*pi/beta) * repmat((1:M),[N_baths,1]) ;

end