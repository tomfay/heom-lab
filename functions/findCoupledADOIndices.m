function [lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = findCoupledADOIndices(ado_indices)
% finds the row indices and column indices corresponding of the coupled
% ADOs int he truncated set of ADOs, specified such that the n_jk values of
% the ADO in the row index is always less than or equal to those in the
% column index

n_ados = size(ado_indices,1) ;
n_modes = size(ado_indices,2) ;
lower_indices = [] ;
upper_indices = [] ;
coupled_mode_indices = [] ;
truncated_coupled_modes = zeros([n_ados,n_modes],'logical') ;

for J = 1:n_ados
   for k = 1:n_modes
      % for each ado with indices specified in ado_indices, generate the
      % ado indices coupled above it in the heirarchy, and if it exists in
      % the list, add that coupling
      ado_indices_k_plus = ado_indices(J,:) ;
      ado_indices_k_plus(k) = ado_indices_k_plus(k)+1 ;
      [is_in_ado_set,K] = ismember(ado_indices_k_plus,ado_indices((J+1):end,:),'rows') ;
      if (is_in_ado_set) 
         lower_indices = [lower_indices; J] ;
         upper_indices = [upper_indices; K+J] ;
         coupled_mode_indices = [coupled_mode_indices; k] ;
      else 
          truncated_coupled_modes(J,k) = true ;
      end
   end
end


end