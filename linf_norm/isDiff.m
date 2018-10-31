function bool = isDiff(A,x)
% returns true if function (with infinity norm) is differentiable at x

[n,m] = size(A) ;

bool = true ;
for i=1:m
    z = abs(x-A(:,i)) ;
    max = z(1) ;
    bool_i = true ;
    for j=2:n
        if z(j) > max
            max = z(j) ;
            bool_i = true ;
        elseif z(j) == max
            bool_i = false ;
        end
    end
    if bool_i == false
        bool = false ;
        return
    end
end
            
            
            
            