function dny_dxn = calculateDerivative(y,dx,deriv_order,accuracy)

d = size(y,1) ;
n_x = size(y,2) ;
dny_dxn = zeros([d,n_x]) ;

if deriv_order==1
    % second order approximation to first derivative
    if accuracy == 2
        % derivative at start and end has to be calculated using
        % forward/backward approximations
        D_f = (1/dx)*[-1.5,2,-0.5] ;
        D_b = (1/dx)*[0.5,-2,1.5] ;
        D_c = (1/dx)*[-0.5,0,0.5] ;
        dny_dxn(:,1) = sum(y(:,1:3).*D_f,2) ;
        for k = 2:(n_x-1)
           dny_dxn(:,k) = sum(y(:,(k-1):(k+1)).*D_c,2) ; 
        end
        dny_dxn(:,n_x) = sum(y(:,(n_x-2):n_x).*D_b,2) ;
    % fourth order approximation to first derivative
    elseif accuracy == 4     
        % derivative at start and end has to be calculated using
        % forward/backward approximations
        D_f = (1/dx)*[-25/12, 4, -3, 4/3, -1/4] ;
        D_b = - D_f(5:-1:1) ;
        D_c = (1/dx)*[1/12,-2/3,0,2/3,-1/12] ;
        dny_dxn(:,1) = sum(y(:,1:5).*D_f,2) ;
        dny_dxn(:,2) = sum(y(:,2:6).*D_f,2) ;
        for k = 3:(n_x-2)
           dny_dxn(:,k) = sum(y(:,(k-2):(k+2)).*D_c,2) ; 
        end
        dny_dxn(:,n_x-1) = sum(y(:,(n_x-5):(n_x-1)).*D_b,2) ;
        dny_dxn(:,n_x) = sum(y(:,(n_x-4):n_x).*D_b,2) ;
    end
    
elseif deriv_order == 2
    % second derivative at second order accuracy
    if accuracy == 2
        % derivative at start and end has to be calculated using
        % forward/backward approximations
        D_f = (1/dx)^2*[2,-5,4,-1] ;
        D_b = (1/dx)^2*[-1,4,-5,2] ;
        D_c = (1/dx)^2*[1,-2,1] ;
        dny_dxn(:,1) = sum(y(:,1:4).*D_f,2) ;
        for k = 2:(n_x-1)
           dny_dxn(:,k) = sum(y(:,(k-1):(k+1)).*D_c,2) ; 
        end
        dny_dxn(:,n_x) = sum(y(:,(n_x-3):n_x).*D_b,2) ;
        
    % fourth order approximation to second derivative
    elseif accuracy == 4     
        % derivative at start and end has to be calculated using
        % forward/backward approximations
        D_f = (1/dx)^2*[15/4,-77/6,107/6,-13,61/12,-5/6] ;
        D_b = D_f(6:-1:1) ;
        D_c = (1/dx)^2*[-1/12,4/3,-5/2,4/3,-1/12] ;
        dny_dxn(:,1) = sum(y(:,1:6).*D_f,2) ;
        dny_dxn(:,2) = sum(y(:,2:7).*D_f,2) ;
        for k = 3:(n_x-2)
           dny_dxn(:,k) = sum(y(:,(k-2):(k+2)).*D_c,2) ; 
        end
        dny_dxn(:,n_x-1) = sum(y(:,(n_x-6):(n_x-1)).*D_b,2) ;
        dny_dxn(:,n_x) = sum(y(:,(n_x-5):n_x).*D_b,2) ;
    end  
end

end