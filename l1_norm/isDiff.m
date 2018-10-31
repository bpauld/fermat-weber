function bool = isDiff(A,x)
%param A : nxm matrix containing anchor points
%param x : nx1 vector
%returns true if the function is differentiable at x, false
%otherwise ( function with l1 norm)

[n,m] = size(A) ;

bool = true ;
for i=1:m
    for j=1:n
        if x(j) == A(j,i)
            bool = false ;
            return
        end
    end
end