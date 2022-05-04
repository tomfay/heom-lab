function [ado_indices,ado_gammas,lower_indices,upper_indices,coupled_mode_indices,truncated_coupled_modes,ado_indices_term,modes_term,term_indices] = generateHierarchyDepthTrunc(L_max,nus)
% Generates the hierarchy up to level L_max, including all couplings and
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
ado_indices_term = [] ;
modes_term = {} ;
term_indices = {} ;


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
          % if it is a new index set, add it to the list and store the
          % coupling
          if ~is_member
             ado_indices_L = [ado_indices_L; ado_indices_new] ;
             ado_gammas = [ado_gammas;sum(ado_indices_new.*nus)] ;
             n_L = n_L + 1 ;
             lower_indices = [lower_indices; J_r] ;
             K = n_below + n_L ;
             upper_indices = [upper_indices; K] ;
             coupled_mode_indices = [coupled_mode_indices ; k] ;
          % if it isn't new, just store the coupling information
          else
             lower_indices = [lower_indices; J_r] ;
             K = index + n_below ;
             upper_indices = [upper_indices; K] ;
             coupled_mode_indices = [coupled_mode_indices ; k] ;
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
truncated_coupled_modes = zeros([n_ados,n_modes],'logical') ;
truncated_coupled_modes((end-n_L+1):end,:) = true ;

% construct the ADO indices for the terminating layer
% only ados in level L_max get truncated, in all modes, so store this
% information here
truncated_coupled_modes((end-n_L+1):end,:) = true ;
n_below = size(ado_indices,1) ;
n_lower = size(ado_indices_lower,1) ;
ado_indices_L = zeros([0,n_modes]) ;
ado_indices_L_term = zeros([0,n_modes]) ;
n_below_term = size(ado_indices_term,1) ;
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
        gamma_ado = sum(ado_indices_new.*nus) ;
        if (~is_member_term )
            ado_indices_L_term = [ado_indices_L_term; ado_indices_new] ;
            J_term = size(ado_indices_L_term,1) + n_below_term ;
            modes_term{J_term} = [k] ;
            term_indices{J_term} = [J_r] ;
        elseif (is_member_term)
            J_term = index_term + n_below_term ;
            modes_term{J_term} = [modes_term{J_term},k] ;
            term_indices{J_term} = [term_indices{J_term},J_r] ;
        end
    end
end
ado_indices = [ado_indices;ado_indices_L] ;
ado_indices_term = [ado_indices_term;ado_indices_L_term] ;

end