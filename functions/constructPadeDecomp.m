function [nus,cs,cbars,Delta] = constructPadeDecomp(omega_D,lambda_D,beta,N,approximant_type)
if (approximant_type == "[N-1/N]")
    ns = (1:(2*N))' ;
    Lambda = full(spdiags([1./sqrt((2*ns+1).*(2*(ns+1)+1)),1./sqrt((2*ns+1).*(2*(ns-1)+1))],[-1,1],2*N,2*N)) ;
    ns = (1:(2*N-1))' ;
    Lambda_tilde = full(spdiags([1./sqrt((2*ns+3).*(2*(ns+1)+3)),1./sqrt((2*ns+3).*(2*(ns-1)+3))],[-1,1],2*N-1,2*N-1)) ;

    eigvals_Lambda = eig(Lambda) ;
    eigvals_Lambda_tilde = eig(Lambda_tilde) ;
    zeta = transpose(2./eigvals_Lambda_tilde((N+1):end)) ;
    zeta = sort(zeta) ;
    zeta_sq = zeta.*zeta ;
    xi = transpose(2./eigvals_Lambda((N+1):end)) ;
    xi = sort(xi) ;
    xi_sq = xi.*xi ;
    eta = zeros([1,N]) ;
    for j = 1:N
        num_inds = [1:(N-1)] ;
        denom_inds = [1:(j-1),(j+1):(N)] ;
        num = zeta_sq(num_inds) - xi_sq(j) ;
        denom = xi_sq(denom_inds) - xi_sq(j) ;
        eta(j) = (N^2 +(1.5*N)) * prod(num)/prod(denom) ;
    end
elseif (approximant_type == "[N/N]")
    ns = (1:(2*N+1))' ;
    Lambda = full(spdiags([1./sqrt((2*ns+1).*(2*(ns+1)+1)),1./sqrt((2*ns+1).*(2*(ns-1)+1))],[-1,1],2*N+1,2*N+1)) ;
    ns = (1:(2*N))' ;
    Lambda_tilde = full(spdiags([1./sqrt((2*ns+3).*(2*(ns+1)+3)),1./sqrt((2*ns+3).*(2*(ns-1)+3))],[-1,1],2*N,2*N)) ;

    eigvals_Lambda = eig(Lambda) ;
%     eigvals_Lambda
    eigvals_Lambda_tilde = eig(Lambda_tilde) ;
%     eigvals_Lambda_tilde 
    zeta = transpose(2./eigvals_Lambda_tilde((N+1):end)) ;
    zeta = sort(zeta) ;
    zeta_sq = zeta.*zeta ;
    xi = transpose(2./eigvals_Lambda((N+2):end)) ;
    xi = sort(xi) ;
    xi_sq = xi.*xi ;
    eta = zeros([1,N]) ;
    for j = 1:N
        num_inds = [1:(N)] ;
    denom_inds = [1:(j-1),(j+1):(N)] ;
    num = zeta_sq(num_inds) - xi_sq(j) ;
    denom = xi_sq(denom_inds) - xi_sq(j) ;
        R_N = 1/(4*(N+1)*(2*(N+1)+1)) ; 
        eta(j) = (0.5*R_N) * prod(num)/prod(denom) ;
    end
end

nus = [omega_D,xi/beta];
cs = zeros([1,N+1]) ;
cs(1) = (2*lambda_D/beta) * (1 - sum(2*eta.*omega_D*omega_D./((xi/beta).^2 - omega_D^2)) -0.5i * omega_D*beta) ;
cs(2:(N+1)) = (4*lambda_D/beta) .* eta.*omega_D.*(xi/beta)./((xi/beta).^2 - omega_D^2) ;
cbars = conj(cs) ;
if (approximant_type == "[N-1/N]")
    Delta = 0 ;
elseif (approximant_type == "[N/N]")
% low temp correction term
    R_N = 1/(4*(N+1)*(2*(N+1)+1)) ;
    Delta = 2*lambda_D*beta*omega_D * R_N ;
end

end