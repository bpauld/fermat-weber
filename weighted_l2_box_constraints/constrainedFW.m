function [output,nbOfIterations,cpu_time] = constrainedFW(H,boundaries,omega,A,eps,kmax,method)
%param H : matrix defining inner product
%param boundaries : n x 2 matrix defining the boundaries of the rectangle
%constraints
%param omega : 1 x m matrix defining the weights
%param A : n x m matrix defining the anchor points

% output = 0 if minimal anchor point is solution
% output = 1 if maximum number of iterations is exceeded
% output = 2 otherwise

[n,m] = size(A) ;
S = sqrtm(H) ;
x0 = 0 ;


%generate anchor points not on the boundary
%{
A = round(100*rand(n,m));
%replace anchor points on the boundary
for i=1:m
    for j=1:n
        if (A(j,i) == boundaries(j,1) || A(j,i) == boundaries(j,2))
            onBound = true ;
            while onBound
                A(j,i) = 100*rand ;
                if (A(j,i) ~= boundaries(j,1) && A(j,i) ~= boundaries(j,2))
                    onBound = false ;
                end
            end
        end
    end
end
%}



%find the anchor points with the minimum fucntion value within and outside
%the boundaries

fp = 0 ;
for i=1:m
    fp = fp + omega(i)*norm(S*(A(:,1)-A(:,i))) ;
end

p=1;
fmin = fp ;

j = -1 ;
fj = -1 ;

if withinBoundaries(A(:,1),boundaries)
    j = 1 ;
    fj = fp ;
end

for i = 2:m
   fp = 0;
   for k = 1:m
      fp = fp + omega(k)*norm(S*(A(:,k)-A(:,i)));
   end
   if fp < fmin 
       p = i ;
       fmin = fp ;
   end
   if (withinBoundaries(A(:,i),boundaries) && fp < fj) || (withinBoundaries(A(:,i),boundaries) && fj == -1)
       j = i ;
       fj = fp ;
   end
end





if j==p
    %in that case the anchor with minimum function value is in the
    %rectangle, so we check the condition from propostion 2.2
    crit = zeros(n,1) ;
    for i=1:m
        if i~=p
            num = S*(A(:,p) - A(:,i)) ;
            crit = crit + omega(i)*num/norm(num) ;
        end
    end
    normCrit = norm(crit) ;
    
    if normCrit <= omega(p)
        %disp('Minimal anchor point is solution') ;
        %A(:,p)
        output = 0 ;
        nbOfIterations = 0 ;
        cpu_time = 0 ;
        return
    end
end
if j==-1
    %set x0 to be a random point within the rectangle
    x0 = zeros(n,1) ;
    for i=1:n
        x0(i) = (boundaries(i,1)+boundaries(i,2))/2 ;
    end
    if withinBoundaries(x0,boundaries) == false
        disp('problem with x0') ;
    end
else
    %   Compute descent direction
    Rj = zeros(n,1);
    for i = 1:m
        if i ~= j
            adiff = A(:,j) - A(:,i);   
            Rj = Rj + omega(i)*adiff/norm(S*adiff);
        end
    end
    normRj = norm(S*Rj);

    dj = -Rj / normRj ;
        
    t = 2 ;
    x = A(:,j) ;
    xnew = x + t*dj ;
    fxnew = 0 ;
    for i =1:m
        fxnew = fxnew + omega(i)*norm(S*(xnew-A(:,i))) ;
    end
    while (fxnew >= fj || ~(withinBoundaries(xnew,boundaries)))
        t= t/2 ;
        
        xnew = x + t*dj ;
        fxnew = 0 ;
        for i =1:m
            fxnew = fxnew + omega(i)*norm(S*(xnew-A(:,i))) ;
        end
    end
    x0 = xnew ;   
end
fx0 = 0 ;


if method == 'newton'
    tic
    [x,nbOfIterations] = projectedNewton(x0,H,boundaries,A,omega,kmax) ;
    cpu_time = toc ;
end
if method == 'fista'
    tic
    [x,nbOfIterations] = fista(A,omega,H,x0,boundaries,eps,kmax) ;
    cpu_time = toc ;
end

if nbOfIterations >= kmax
    output = 1 ;
else
    output = 2 ;    
end
end



                
