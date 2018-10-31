function values = subgradient_values(A,omega,p,nbOfValues,boundaries)
%param A : n x m matrix. Each of the m columns represent an anchor point.
%param omega : 1 x m vector representing weight of each anchor point
%param boundaries (optionnal) : n x 2 matrix representing the boundaries of the set
%if boundaries is not given, then the method is done on R^n
%param p : norm for the FW problem


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

values = zeros(nbOfValues,1) ;


%compute function value at xk
fk = FW(A,omega,xk,p) ;

values(1,1) = fk ;

%f_best is the smallest function value attained so far
f_best = fk ;
x_best = x0 ;
k_best = 0 ;

for k=2:nbOfValues
    % compute one subgradient at xk
    if p==inf
        sub_fk = linf_subgradient(A,omega,xk) ;
    elseif p==1
        sub_fk = l1_subgradient(A,omega,xk) ;
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
    fx_new = FW(A,omega,x_new,p) ;
    
    values(k,1) = fx_new ;
    
    %update f_best if necessary
    if fx_new < f_best
        f_best = fx_new ;
        x_best = x_new ;
        k_best = k+1 ;
    end
    
    
    
    xk = x_new ;
    fk = fx_new ;
    
end