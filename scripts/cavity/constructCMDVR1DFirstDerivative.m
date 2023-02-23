function D1 = constructCMDVR1DFirstDerivative(dx,n_DVR)
    % constructs the Colbert-Miller DVR representation of the 1D d^2/dx^2
    % operator with n_DVR DVR points and a spacing between points of dx

    % indices of the DVR points
    inds = 1:n_DVR ;
    % off-diagonal elements (-2/dx^2) (-1)^(n-m)/((n-m)^2)
    D1 = -(1/(dx))*((-1).^(inds-inds'))./((inds-inds'));
    % diagonal elements (-pi^2/(3 dx^2))
    D1(logical(eye(n_DVR))) = 0 ;

end