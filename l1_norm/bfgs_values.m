function function_value = bfgs_values(A,omega,x0,c1,c2,nbOfValues)
%implementation of bfgs method with inexact line search for nonsmooth
%optimization with starting matrix H0 = I
%param A : n x m matrix. Each of the m columns represent an anchor point.
%param omega : 1 x m vector representing weight of each anchor point
%param x0: nx1 vector. Starting point of the algorithm. The underlying
%function should be differentiable at x0
%params c1 and c2 : params for inexact line search, should be such that
%0 < c1 < c2 < 1
%nbOfValues = number of iterations wanted

%output function_value : vector of size nbOfValues containing function
%values at each iteration

[n,m] = size(A) ;
H = eye(n) ;

function_value = zeros(nbOfValues,1) ;

xk = x0 ;

gk = l1_gradient(A,omega,xk) ;
fxk = FW(A,omega,xk,1) ;

function_value(1,1) = fxk ;

for k=2:nbOfValues
    
    d =  - H\gk ;
    
    t = inexactLineSearch(xk,gk,A,omega,d,c1,c2) ;
    
    %update
    x_new = xk + t*d ;
    sk = x_new - xk ;
    
    %compute function value at x_new
    fx_new = FW(A,omega,xk,1) ;
    
    %keep value 
    function_value(k,1) = fx_new ;
    
    %compute gradient at x_new
    g_new = l1_gradient(A,omega,x_new) ;
    
    yk = g_new - gk ;
    H = H + yk*yk'/(yk'*sk) - H*(sk*sk')*H/(sk'*H*sk) ;
    
    gk = g_new ;
    xk = x_new ;
    fxk = fx_new ;
end
end


