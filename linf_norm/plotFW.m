function plotFW(A,omega,p,x,y)
[t,k1] = size(x) ;
[t,k2] = size(y) ;
Z = zeros(k1,k2) ;
for i =1:k1
    for j=1:k2
        Z(i,j) = FW(A,omega,[x(i);y(j)],p) ;
    end
end

figure
[X,Y] = meshgrid(x,y) ;
[C,h] = contour(X,Y,Z,30);   
clabel(C,h)