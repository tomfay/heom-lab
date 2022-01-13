function M_max = findMaxMjWeightedCutoffBO(L_cut, Omega, lambda, gamma, beta)
% finds the maximum values of M_j, max number of matsurbara frequencies for
N_baths = length(gamma_D) ;
M_baths = zeros([1,N_baths]) ;
for j = 1:N_baths
    a = 2*pi/(beta*gamma(j)) ;
    b = beta*Omega(j)/(2*pi) ;
    c = L_cut *sqrt(lambda(j)*Omega(j)*Omega(j)/(8*pi*gamma(j))) * (beta/(2*pi)) ;
    M_guess = 0 ;
    f = @(x) sqrt(x.*abs( x.*x - a*(x.*x+ b ).^2 ))- c ;
    M_baths(j) = floor(fzero(f,M_guess)) ;
end
M_max = max(M_baths) ;

end