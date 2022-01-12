function [ado_indices,ado_gammas] = generateADOIndices(Gamma_cut,nus)
% generate a list of n_00,n_01,... n_NbM excitation indices thata re
% included in the hierarchy

% nus is an array of nu_jk (dimensions N_baths x (M+1))
N_baths = size(nus, 1) ; % number of baths
M = size(nus,2)-1 ; % number of matsubara modes included in the expansion
N_inds = N_baths*(M+1) ;

% first set to be included is the set where all n_jk = 0
ado_indices = zeros([1,N_inds]) ;
ado_gammas = [0] ;


% find max n_jk for each mode,
n_maxs = floor(Gamma_cut./nus) ;
n_maxs = reshape(n_maxs',[1,N_inds]) ;
N_HEOM_max = prod(n_maxs+1,'all') ;
nus_row = reshape(transpose(nus),[1,N_inds]) ;

% (n_maxs+1).*nus_row
% simple - probably not very efficient implementation to generate all sets
% of ado indices to include in the heirachy
for J = 2:N_HEOM_max
    n_jks = getTensorFromVectorIndex(J,n_maxs+1)-1 ;
    gamma_n_jk = sum(nus_row.*n_jks) ;
    if gamma_n_jk < Gamma_cut 
       ado_indices = [ado_indices; n_jks]  ;
       ado_gammas = [ado_gammas ; gamma_n_jk] ;
    end
    
end

% % a more efficient algorithm for generating the ADO list, along with
% % connectivities
% n_checked = 0 ;
% n_ados = size(ado_indices,1) ;
% while 
% 
% end