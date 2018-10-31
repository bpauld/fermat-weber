function t = withinBoundaries(vector,boundaries)
[n,k] = size(boundaries) ;
%k should equal 2
t = false ;

for i=1:n
    if vector(i) < boundaries(i,1) || vector(i) > boundaries(i,2)
        return
    end
end
t = true ;
end