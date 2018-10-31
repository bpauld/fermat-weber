function A = createAnchorPointsNotOnBound(m,boundaries)
[n,k] = size(boundaries) ;
A = 100*rand(n,m) ;
for j =1:m
    for i=1:n
        while A(i,j) == boundaries(i,1) || A(i,j) == boundaries(i,2)
            A(i,j) = 100*rand ;
        end
    end
end
end