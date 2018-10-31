function sub = linf_subgradient(A,omega,x)
% computes the subgradient of the FW function with l-infinity norm at x

[n,m] = size(A) ;
sub = zeros(n,1) ;

for i=1:m
    sub_i = zeros(n,1) ;
    [maximum,argmax] = max (abs(x-A(:,i))) ;
    if x(argmax) > A(argmax,i)
        sub_i(argmax) = omega(i) ;
    elseif x(argmax) < A(argmax,i)
        sub_i(argmax) = -omega(i) ;
    else
        sub_i(argmax) = 0 ;
    end
    sub = sub + sub_i ;
end
end