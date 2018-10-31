m = 100 ;
m

H = eye(2) ;
for i=-4:4
    H(2,2) = 10^(i) ;
    generalizedFW(H,100,m) ;
end

H = eye(3) ;
for i=-4:4
    H(3,3) = 10^(i) ;
    generalizedFW(H,100,m) ;
end

H = eye(4) ;
for i=-4:4
    H(4,4) = 10^(i) ;
    generalizedFW(H,100,m) ;
end


H = eye(6) ;
for i=-4:4
    H(6,6) = 10^(i) ;
    generalizedFW(H,100,m) ;
end

H = eye(8) ;
for i=-4:4
    H(8,8) = 10^(i) ;
    generalizedFW(H,100,m) ;
end

H = eye(10) ;
for i=-4:4
    H(10,10) = 10^(i) ;
    generalizedFW(H,100,m) ;
end