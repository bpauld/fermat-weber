%clear

m = 10 ;

for n = 3
    H = eye(n) ;
    eps = [0.0001 ; 0.001 ; 0.01 ; 0.1 ;1; 10; 100; 1000; 10000] ;
    for t = 1:size(eps)
        H(n,n) = eps(t) ;
        avIter = 0 ;
        nbOf0 = 0 ;
        nbOf1 = 0 ;
        nbOf2 = 0 ;
        av_cpu_time = 0 ;
        for k=1:100
            boundaries = 100*rand(n,2) ;
            boundaries = rearrangeBoundaries(boundaries) ;
            A = createAnchorPointsNotOnBound(m,boundaries) ;
            omega = 100*rand(1,m) ;
            [output,nbOfIterations, cpu_time] = constrainedFW(H,boundaries,omega,A) ;
            if output == 0
                nbOf0 = nbOf0 + 1 ;
            elseif output == 1
                nbOf1 = nbOf1 + 1;
            elseif output == 2
                nbOf2 = nbOf2 + 1;
                avIter = avIter + nbOfIterations ;
                av_cpu_time = av_cpu_time + cpu_time ; 
            end
        end
        avIter = avIter/nbOf2 ;
        av_cpu_time = av_cpu_time/nbOf2 ;
            
        fprintf('for n = %4.0f and m = %4.0f \neps = %16.8f\nNumber of Success= %4.0f\nNumber of failures=%4.0f\nAverage number of iterations = %16.4f\nAverage cpu time = %16.4f\n\n',n,m,eps(t),nbOf2,nbOf1,avIter,av_cpu_time);
    end
end
