function values = bfgs_inf_values(A,omega,x0,c1,c2,nbOfValues)
%implementation of bfgs method with inexact line search for nonsmooth
%optimization with starting matrix H0 = I
%param A : n x m matrix. Each of the m columns represent an anchor point.
%param omega : 1 x m vector representing weight of each anchor point
%param x0: nx1 vector. Starting point of the algorithm. The underlying
%function should be differentiable at x0
%params c1 and c2 : params for inexact line search, should be such that
%0 < c1 < c2 < 1

[n,m] = size(A) ;
H = eye(n) ;

values = zeros(nbOfValues,1) ;


xk = x0 ;

%compute gradient at x0
gk = linf_gradient(A,omega,xk) ;

%compute function value
fxk = 0 ;
for i=1:m
    fxk = fxk + omega(i)*norm(xk-A(:,i),inf) ;
end

values(1,1) = fxk ;

for k=2:nbOfValues
        
    
    
    
    d =  - H\gk ;
    
    t = inexactLineSearch(xk,gk,A,omega,d,c1,c2) ;
    
    %update
    x_new = xk + t*d ;
    sk = x_new - xk ;
    
    %compute function value at x_new
    fx_new = 0 ;
    for i=1:m
        fx_new = fx_new + omega(i)*norm(x_new - A(:,i),inf) ;
    end
    
    values(k,1) = fx_new ;
    
    
    %compute gradient at x_new
    g_new = linf_gradient(A,omega,x_new) ;
    
    yk = g_new - gk ;
    H = H + yk*yk'/(yk'*sk) - H*(sk*sk')*H/(sk'*H*sk) ;
    
    gk = g_new ;
    xk = x_new ;
end