function grad = linf_gradient(A,omega,x)
%computes gradient of f (infinity norm) at x, assuming that f is
%differentiable at x

[n,m] = size(A) ;
grad = zeros(n,1) ;

for i = 1:m 
    [maximum,argmax] = max(abs(x-A(:,i))) ;
    grad_i = zeros(n,1) ;
    grad_i(argmax) = omega(i)*sign(x(argmax)-A(argmax,i)) ;
    grad = grad + grad_i ;
end
end
    