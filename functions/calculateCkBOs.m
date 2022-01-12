function c_BO = calculateCkBOs(gamma,Omega,beta,lambda,k)
% calculates the matsubara term coefficients for the brownian oscillator
nu = (2*pi/beta)*k  ;
c_BO = (lambda * gamma * Omega * Omega/(4.0*beta))./(gamma*gamma*nu.*nu - (nu.*nu + Omega * Omega).^2) ; 

end