function [x,numberOfIterations,x_best,k_best,stopping_cond] = projectedSubgradient(A,omega,eps1,eps2,boundaries)
%param A : n x m matrix. Each of the m columns represent an anchor point.
%param omega : 1 x m vector representing weight of each anchor point
%param eps : stopping criterion, algorithm stops when f(x(k+1)) - f(xk) <
%eps
%param boundaries (optionnal) : n x 2 matrix representing the boundaries of the set
% if boundaries is not given, then the method is done on R^n

[n,m] = size(A) ;

bound = true ;
if ~exist('boundaries','var')
    bound = false ;
end


%pick random starting point x0 within the boundaries

x0 = zeros(n,1) ;
if bound
    for i=1:n
        x0(i) = (boundaries(i,1) + boundaries(i,2)) / 2 ;
    end
end



%prepare for loop
k = 0 ;
xk = x0 ;
cond = true ;


%compute function value at xk
fk = 0 ;
for i=1:m
    fk = fk + omega(i)*norm(xk-A(:,i),1) ;
end


%f_best is the smallest function value attained so far
f_best = fk ;
x_best = x0 ;
k_best = 0 ;


while cond
    % compute one subgradient at xk
    sub_fk = zeros(n,1) ;
    for i=1:m
        gi = zeros(n,1) ;
        for j=1:n
            if xk(j) > A(j,i)
                gi(j) = omega(i) ;
            elseif xk(j) < A(j,i)
                gi(j) = - omega(i) ;
            else
                gi(j) = 0 ;
            end
        end
        sub_fk = sub_fk + gi ;
    end
    
    if norm(sub_fk) < eps1
        x = xk ;
        numberOfIterations = k ;
        stopping_cond = 1 ;
        return
    end
    

    
    %compute stepsize
    tk = 1/(sqrt(k+1)) ;
    
    %compute new vector
    x_new = zeros(n,1) ;
    if bound
        x_new = projection(xk-tk*sub_fk,boundaries) ;
    else
        x_new = xk - tk*sub_fk ;
    end
    
    %compute function value at new vector
    fx_new = 0 ;
    for i=1:m
        fx_new = fx_new + omega(i)*norm(x_new-A(:,i),1) ;
    end
    
    %update f_best if necessary
    if fx_new < f_best
        f_best = fx_new ;
        x_best = x_new ;
        k_best = k+1 ;
    end
    
    
    
    if abs(fx_new - fk) < eps2
        x = x_new ;
        numberOfIterations = k + 1 ;
        stopping_cond = 2 ;
        cond = false ;
    else
        k = k+1 ;
        xk = x_new ;
        fk = fx_new ;
    end
    
    
end

