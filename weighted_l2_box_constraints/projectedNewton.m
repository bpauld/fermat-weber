function [xfinal,nbOfSteps,iterates] = projectedNewton(x0,H,boundaries,A,omega,kmax)
[n,m] = size(A) ;

S = sqrtm(H) ;
k = 0 ;
eps = 0.0001 ;
xk = x0 ;
beta = 0.5 ;
sigma = 0.0000025 ;


iterates = zeros(n,2001) ;



cond = true ;
while cond
    
    
    
    grad = zeros(n,1) ;
    for i=1:m
        grad = grad + omega(i)*H*(xk-A(:,i))/norm(S*(xk-A(:,i))) ;
    end
    
    wk = norm(xk-projection(xk-grad,boundaries)) ;
    eps_k = min(eps,wk) ;
    
    
    
    %Ik is a nx1 vector similar to equation (32). If i is in Ik,
    % then Ik(i) = 1, otherwise Ik(i) = 0
    Ik = zeros(n,1) ;
    for i=1:n
        if ((boundaries(i,1) + eps_k >= xk(i) && grad(i) > 0) || (boundaries(i,2) - eps_k <= xk(i) && grad(i) < 0))
            Ik(i) = 1 ;
        end
    end
    
    hessian = zeros(n,n) ;
    for i=1:m
        vecDiff = xk - A(:,i) ;
        normVecDiff = norm(S*vecDiff) ;
        hessian = hessian + omega(i)*(1/normVecDiff^3)*(normVecDiff^2*H - H*(vecDiff*vecDiff')*H) ;
    end
    Hk = zeros(n,n) ;
    for i=1:n
        for j=1:n
            if i~=j && (Ik(i) == 1 || Ik(j) == 1)
                Hk(i,j) = 0 ;
            else
                Hk(i,j) = hessian(i,j) ;
            end
        end
    end
    %Dk = inv(Hk) ;
    
    
    
    
    pk = Hk\grad ;
    fxk = 0 ;
    mk = 0 ;
    alpha_k = beta^mk ;
    for i=1:m
        fxk = fxk + omega(i)*norm(S*(xk-A(:,i))) ;
    end
    
    xk_alpha = projection(xk-alpha_k*pk,boundaries) ;
    
    fxk_alpha = 0 ;
    for i=1:m
        fxk_alpha = fxk_alpha + omega(i)*norm(S*(xk_alpha-A(:,i))) ;
    end
    
    
    sum_active = 0 ;
    sum_free = 0 ;
    for i=i:n
        if Ik(i) == 0
            sum_free = sum_free + grad(i)*pk(i) ;
        else
            sum_active = sum_active + grad(i)*(xk(i) - xk_alpha(i)) ;
        end
    end
    
    while fxk - fxk_alpha  < sigma * ( alpha_k * sum_free + sum_active)
        mk = mk + 1 ;
        alpha_k = beta^mk ;
        xk_alpha = projection(xk-alpha_k*pk,boundaries) ;
    
        fxk_alpha = 0 ;
        for i=1:m
            fxk_alpha = fxk_alpha + omega(i)*norm(S*(xk_alpha-A(:,i))) ;
        end
        
        sum_active = 0 ;
        for i=i:n
            if Ik(i) == 1
                sum_active = sum_active + grad(i)*(xk(i) - xk_alpha(i)) ;
            end
        end
    end
    
    if (norm(xk-xk_alpha) < 1e-8 ) || k > kmax
        %if (mk == 0 && norm(xk - xk_alpha) < eps_prime) || k > kmax
        %in that case we are done
        cond = false ;
    end
    

    
    
    xk = xk_alpha ;
    k = k+1 ;
    

    iterates(:,k) = xk ;
    
end   
xfinal = xk ;
nbOfSteps = k ;
end