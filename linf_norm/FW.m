function f = FW(A,omega,x,p)
%function that computes the Fermat Weber function value at x
%param A : n x m matrix. Each of the m columns represent an anchor point.
%param omega : 1 x m vector representing weight of each anchor point
%param p : norm wanted, must format to norm function provided by Matlab

[n,m] = size(A) ;

f = 0 ;
for i=1:m
    f = f + omega(i)*norm(x-A(:,i),p) ;
end
end
