function [ado_indices,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateADOIndexSet(L_max,N_modes)

ado_indices = zeros([1,N_modes]) ;
ado_indices_lower = ado_indices ;
for L = 1:L_max
   ado_indices_L = zeros([0,N_modes]) ;
   n_lower = size(ado_indices_lower,1) ;
   % generate new sets of excitation numbers from the lower set
   for k = 1:N_modes
       for r = 1:n_lower
          ado_indices_r = ado_indices_lower(r,:) ;
          % generate a new set of indices from a set for lower L
          ado_indices_new = ado_indices_r ;
          ado_indices_new(k) = ado_indices_new(k) + 1 ;
          [is_member,index] = ismember(ado_indices_new,ado_indices_L,'rows') ;
          if ~is_member
             ado_indices_L = [ado_indices_L; ado_indices_new] ;
          end
       end
   end
   % store the new set
   ado_indices = [ado_indices;ado_indices_L] ;
   ado_indices_lower = ado_indices_L ;
end

end