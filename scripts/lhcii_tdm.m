% a script for obtaining the LHCII transition dipole moment operators
chl_data = readlines('./scripts/data/1rwt-chainc-chl.pdb') ;
r_NB = [] ;
r_ND = [] ;
inds_chla = [] ;
inds_chlb = [] ;
n_lines = size(chl_data,1) ;
n_NB = 0 ;
n_ND = 0 ;

for k = 1:n_lines
    strings = split(chl_data(k)) ;
    if (strings(2)=="NB")
        n_NB = n_NB + 1 ;
        r = [str2double(strings(6)),str2double(strings(7)),str2double(strings(8))] ;
        r_NB(n_NB,:) = r ;
        if (strings(3)=="CHL")
            inds_chlb = [inds_chlb,n_NB] ;
        end
        if (strings(3)=="CLA")
            inds_chla = [inds_chla,n_NB] ;
        end
    end
    if (strings(2)=="ND")
        n_ND = n_ND + 1 ;
        r = [str2double(strings(6)),str2double(strings(7)),str2double(strings(8))] ;
        r_ND(n_ND,:) = r ;
    end
end

unit_vectors = r_ND - r_NB ;
unit_vectors = unit_vectors ./ sqrt(sum(unit_vectors.^2,2)) ;

mu_chla = 4.0 ;
mu_chlb = 3.4 ;
mu = unit_vectors ;
mu(chla_inds,:) = mu_chla * mu(chla_inds,:) ;
mu(chlb_inds,:) = mu_chlb * mu(chlb_inds,:) ;

% construct the transition operators (in units of seconds)
Cm_per_D = 3.33564e-30 ;
epsilon_0 = 8.8541878128e-12 ;
c_0 = 299792458 ;
hbar = 1.054571817e-34 ;

% this gives the transition operator pre-factor in s
T = Cm_per_D / sqrt( 6* (pi^2)*(c_0^3)*epsilon_0*hbar  ) ;
omega_0 = 2*pi*c_0 * 15287e2 ;
% Convert T to units of cm D^-1
T = (2*pi * c_0 * 1e2 ) * T ;

% sanity check for T
Gamma_in_wn = 2*T^2 *(pi * (15287)^3)*(mu_chla^2) ;
Gamma_in_Hz = Gamma_in_wn*(2*pi*c_0*1e2) ;

T_op = T * mu ;
