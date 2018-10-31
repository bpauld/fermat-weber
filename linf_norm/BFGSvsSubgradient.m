clear


eps1 = 1e-8 ;
eps2 = 1e-4 ;
c1 = 0.25 ;
c2 = 0.5 ;

m = 10 ;
nbOfTest = 100 ;

for n= 2:2:10
    n
    for t = 1:1
        avIter_bfgs = 0 ;
        counter_bfgs = 0 ;
        avIter_sub = 0 ;
        avIterBest_sub = 0 ;
        counter_sub = 0 ;
        for k=1:nbOfTest
            A = 100*rand(n,m) ;
            omega = 100*rand(1,m) ;
            x0 = zeros(n,1) ;
            [x,nbIt,stopping_cond] = bfgs_inf(A,omega,x0,c1,c2,eps1,eps2) ;
            avIter_bfgs = nbIt + avIter_bfgs ;
            if stopping_cond == 1
                counter_bfgs = counter_bfgs + 1 ;
            end
            [x,nbIt,x_best,k_best,stopping_cond] = subgradientMethod(A,omega,eps1,eps2,inf) ;
            avIter_sub = avIter_sub + nbIt ;
            avIterBest_sub = avIterBest_sub + k_best ;
            if stopping_cond == 1
                counter_sub = counter_sub + 1 ;
            end
        end
        avIter_bfgs = avIter_bfgs/nbOfTest ;
        avIter_sub = avIter_sub/nbOfTest ;
        avIterBest_sub = avIterBest_sub/nbOfTest ;
        avIter_bfgs
        counter_bfgs
        avIter_sub
        avIterBest_sub
        counter_sub
    end
end
            