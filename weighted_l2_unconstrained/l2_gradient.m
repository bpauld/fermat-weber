function grad = l2_gradient(H,A,omega,x)
%computes gradient at x of function defined by A, omega, and with norm induced
%by H.

[n,m] = size(A) ;

S=sqrt(H) ;
grad = zeros(n,1) ;
for i=1:m
    grad = grad + omega(i)*H*(x-A(:,i))/norm(S*(x-A(:,i))) ;
end