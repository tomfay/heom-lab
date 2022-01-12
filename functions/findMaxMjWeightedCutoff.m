function M_max = findMaxMjWeightedCutoff(L_cut, gamma_D, lambda_D, beta)
% finds the maximum values of M_j, max number of matsurbara frequencies for
% a given cut-off frequency a la Dijkstra https://doi.org/10.1063/1.4997433
N_baths = length(gamma_D) ;
M_baths = zeros([1,N_baths]) ;
for j = 1:N_baths
    a = gamma_D(j)*beta/(2*pi) ;
    b = (beta / (2*pi))^3 *  4 * lambda_D(j)*gamma_D(j)*(L_cut^2)/beta ;
    M_guess = 0 ;
    f = @(x) x.*abs(x.*x - (a*a)) - b ;
    M_baths(j) = floor(fzero(f,M_guess)) ;
end
M_baths
M_max = max(M_baths) ;

end