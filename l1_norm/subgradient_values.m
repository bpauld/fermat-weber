function function_values = subgradient_values(A,omega,nbOfValues,boundaries)

[n,m] = size(A) ;

function_values = zeros(nbOfValues,1) ;

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

xk = x0 ;


fk = FW(A,omega,xk,1) ;
function_values(1,1) = fk ;

for k=2:nbOfValues
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
    fx_new = FW(A,omega,x_new,1) ;
    function_values(k,1) = fx_new ;
    
    xk = x_new ;
end
    


