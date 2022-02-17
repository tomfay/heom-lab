function [omegas,weights] = discretiseSpecDenGL(N,rho)
    % Define R(x) 
    R = @(x,a) integral(rho,0,x)-a ;
    [x_vals,weights] = generateQuadratureGL(N,0,1) ; 
    weights = weights' ;
    opts = optimset('Display','notify') ;
    for i = 1:N 
       omega_0 = 0.0 ;
       f = @(x) R(x,x_vals(i)) ;
%        omegas(i) = fsolve(f,omega_0) ;
       omegas(i) = fzero(f,omega_0,opts) ;
       omega_0 = omegas(i) ;
    end

end