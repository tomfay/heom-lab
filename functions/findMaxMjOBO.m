function M_baths = findMaxMjOBO(Gamma_cut, beta)
% finds the maximum values of M_j, max number of matsurbara frequencies for
% a given cut-off frequency a la Dijkstra https://doi.org/10.1063/1.4997433
M_baths = max([floor((beta*Gamma_cut/(2*pi)) ),0]) ;

end