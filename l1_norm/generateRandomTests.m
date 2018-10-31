
eps1 = 1e-8 ;
eps2 = 1e-4 ;

dim = [2;3; 4;5; 6; 8;10] ;
m = 10 ;
nbOfTest = 100 ;
avIter = 0 ;
avBestIter = 0 ;

for j =1:7
    
    n = dim(j) ;
    n
    avIter = 0 ;
    avBestIter = 0 ;
    counter = 0 ;
    for i=1:nbOfTest
        
        A = 100*rand(n,m) ;
        omega = 100*rand(1,m) ;
        x0 = -ones(n,1) ;
        [x,k,stopping_cond] = bfgs(A,omega,x0,0.25,0.5,eps1,eps2) ;
        if stopping_cond == 1 
            counter = counter +1 ;
        end
        avIter = avIter + k ;
    end
    avIter = avIter/nbOfTest ;
    avIter
    counter
end
