function [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes] = generateHierarchyWeightedCutoff2Trunc(epsilon_cut,nus,cs)
% Generates the hierarchy using a frequency based cut-off, 
% the decay rates of each ADO in the hierarchy

% get number of modes
n_modes = numel(nus) ;

% generate the first set of ADO indices corresponding to 
ado_indices = zeros([1,n_modes]) ;
ado_indices_lower = ado_indices ;
ado_gammas = [0] ;

% create empty vectors for things
lower_indices = [] ;
upper_indices = [] ;
coupled_mode_indices = [] ;
truncated_coupled_modes = zeros([1,n_modes],'logical') ;

% calculate L_max for the hierarchy
L_max_modes = ceil((1/epsilon_cut)^2*(abs(cs))./(abs(nus).^2)) ;
L_max = max(L_max_modes) ;

% mode weights in calculating 
abs_cs = abs(cs) ;
sqrt_abs_cs = sqrt(abs_cs) ;
L_max_modes



for L = 1:L_max
   n_below = size(ado_indices,1) ;
   n_lower = size(ado_indices_lower,1) ;
   ado_indices_L = zeros([0,n_modes]) ;
   n_L = 0 ;
  
   % generate new sets of excitation numbers from the lower set
   for k = 1:n_modes
       for r = 1:n_lower
          % get a set of ado indices from the level (L-1)
          ado_indices_r = ado_indices_lower(r,:) ;
          J_r = n_below - n_lower + r ;
          % generate a new set of indices from a set for lower L by adding
          % an additional excitation in mode k
          ado_indices_new = ado_indices_r ;
          ado_indices_new(k) = ado_indices_new(k) + 1 ;
          % check to see if this new set of ado indices has already been
          % generated
          [is_member,index] = ismember(ado_indices_new,ado_indices_L,'rows') ;
          gamma_ado = sum(ado_indices_new.*nus) ;
          cutoff_func = sum(sqrt(ado_indices_new).*sqrt_abs_cs)/abs(gamma_ado) ;
          is_included = (cutoff_func>epsilon_cut) ;
          % if it is a new index set, add it to the list and store the
          % coupling
          if ((~is_member) & is_included)
             ado_indices_L = [ado_indices_L; ado_indices_new] ;
             ado_gammas = [ado_gammas;gamma_ado] ;
             truncated_coupled_modes = [truncated_coupled_modes;zeros([1,n_modes],'logical')] ;
             n_L = n_L + 1 ;
             lower_indices = [lower_indices; J_r] ;
             K = n_below + n_L ;
             upper_indices = [upper_indices; K] ;
             coupled_mode_indices = [coupled_mode_indices ; k] ;
          % if it isn't new, just store the coupling information
          elseif is_member
             lower_indices = [lower_indices; J_r] ;
             K = index + n_below ;
             upper_indices = [upper_indices; K] ;
             coupled_mode_indices = [coupled_mode_indices ; k] ;
          elseif (~is_included)
             truncated_coupled_modes(J_r,k) = true ;
          end
       end
   end
   % store the new set
   ado_indices = [ado_indices;ado_indices_L] ;
   ado_indices_lower = ado_indices_L ;
end

% only ados in level L_max get truncated, in all modes, so store this
% information here
n_ados = size(ado_indices,1) ;
truncated_coupled_modes((end-n_L+1):end,:) = true ;

end