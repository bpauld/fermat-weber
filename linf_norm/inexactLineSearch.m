function t = inexactLineSearch(x,g,A,omega,d,c1,c2)

alpha = 0 ;
beta = 'inf' ;
t = 1 ;
[n,m] = size(A) ;

%compute s, lim sup h(t)/t as t --> 0
%it is symply the directional derivative
s = g'*d ;

%compute fx
fx = 0 ;
for i=1:m
    fx = fx + omega(i) * norm(x-A(:,i),inf) ;
end

cond = true ;

while cond
    
    x_cand = x+t*d ;
    %compute h(t)
    h = 0 ;
    for i =1:m
        h = h + omega(i) * norm(x_cand - A(:,i),inf) ;
    end
    h = h - fx ;
    
    %check if A(t) holds
    if h >= c1*s*t
        beta = t ;
    else
        is_diff = isDiff(A,x_cand) ;
        if is_diff
            %compute gradient at x_cand
            g_cand = linf_gradient(A,omega,x_cand) ;
            h_prime = g_cand'*d ;
            
            if h_prime >c2*s
                return
            end
        end
        alpha = t ;
    end
    
    %update t
    if beta == 'inf'
        t = 2*alpha ;
    else
        t = (alpha + beta ) /2 ;
    end
end   
end