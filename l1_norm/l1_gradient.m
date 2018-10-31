function g = l1_gradient(A,omega,x)
%computes gradient at x

[n,m] = size(A) ;
g = zeros(n,1) ;

for i = 1:m
    gi = sign(x-A(:,i)) ;
    g = g + omega(i)*gi ;
end