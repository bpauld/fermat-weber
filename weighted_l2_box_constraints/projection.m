function y = projection(x,boundaries)
[n,w] = size(x) ;
y = zeros(n,1) ;
for i=1:n
    if x(i) < boundaries(i,1)
        y(i) = boundaries(i,1) ;
    elseif x(i) > boundaries(i,2)
        y(i) = boundaries(i,2) ;
    else
        y(i) = x(i) ;
    end
end
end
            
