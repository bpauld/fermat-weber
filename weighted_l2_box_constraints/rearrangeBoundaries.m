function boundaries = rearrangeBoundaries(bd)
[n,k] = size(bd) ;
boundaries = zeros(n,2) ;
for i=1:n
    if bd(i,1) > bd(i,2)
        boundaries(i,1) = bd(i,2) ;
        boundaries(i,2) = bd(i,1) ;
    else
        boundaries(i,1) = bd(i,1) ;
        boundaries(i,2) = bd(i,2) ;
    end
end
end