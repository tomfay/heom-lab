function [D2,x_DVR] = constructCMDVR1DSecondDerivative(dx,n_DVR)
    % constructs the Colbert-Miller DVR representation of the 1D d^2/dx^2
    % operator with n_DVR DVR points and a spacing between points of dx
    
    % x_DVR - can be shifted arbitrarily
    x_DVR = (( 0:(n_DVR-1))*dx) ;
    
    % indices of the DVR points
    inds = double(1:(int32(n_DVR))) ;
    % off-diagonal elements (-2/dx^2) (-1)^(n-m)/((n-m)^2)
    D2 = (-2/(dx*dx))*((-1).^(inds-inds'))./((inds-inds').*(inds-inds'));
    % diagonal elements (-pi^2/(3 dx^2))
    D2(logical(eye(int32(n_DVR)))) = -(pi*pi/(3.0*dx*dx)) ;
   

end