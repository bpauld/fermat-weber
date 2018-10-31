%clear

eps_stop = 1e-8 ;
kmax = 50000 ;
stopping_criterion = 'function_value' ;

m = 10 ;

for n= 2:2:10 
    H = eye(n) ;
    eps = [0.0001 ; 0.001 ; 0.01 ; 0.1 ;1; 10; 100; 1000; 10000] ;
    for t = 1:size(eps)
        H(n,n) = eps(t) ;
        avIter_newton = 0 ;
        nbOf0_newton = 0 ;
        nbOf1_newton = 0 ;
        nbOf2_newton = 0 ;
        av_cpu_time_newton = 0 ;
        avIter_fista = 0 ;
        nbOf0_fista = 0 ;
        nbOf1_fista = 0 ;
        nbOf2_fista = 0 ;
        av_cpu_time_fista = 0 ;
        for k=1:100
            boundaries = 100*rand(n,2) ;
            boundaries = rearrangeBoundaries(boundaries) ;
            A = createAnchorPointsNotOnBound(m,boundaries) ;
            omega = 100*rand(1,m) ;
            
            %calculate with newton
            [output_newton,nbOfIterations_newton, cpu_time_newton] = constrainedFW(H,boundaries,omega,A,eps_stop,kmax,'newton',stopping_criterion) ;
            if output_newton == 0
                nbOf0_newton = nbOf0_newton + 1 ;
            elseif output_newton == 1
                nbOf1_newton = nbOf1_newton + 1;
            elseif output_newton == 2
                nbOf2_newton = nbOf2_newton + 1;
                avIter_newton = avIter_newton + nbOfIterations_newton ;
                av_cpu_time_newton = av_cpu_time_newton + cpu_time_newton ; 
            end
            
            %now with fista
            [output_fista,nbOfIterations_fista, cpu_time_fista] = constrainedFW(H,boundaries,omega,A,eps_stop,kmax,'fista',stopping_criterion) ;
            if output_fista == 0
                nbOf0_fista = nbOf0_fista + 1 ;
            elseif output_fista == 1
                nbOf1_fista = nbOf1_fista + 1;
            elseif output_fista == 2
                nbOf2_fista = nbOf2_fista + 1;
                avIter_fista = avIter_fista + nbOfIterations_fista ;
                av_cpu_time_fista = av_cpu_time_fista + cpu_time_fista ; 
            end
        end
        
        
        
        avIter_newton = avIter_newton/nbOf2_newton ;
        av_cpu_time_newton = av_cpu_time_newton/nbOf2_newton ;
        
        avIter_fista = avIter_fista/nbOf2_fista ;
        av_cpu_time_fista = av_cpu_time_fista/nbOf2_fista ;
            
        fprintf('Newton :\nfor n = %4.0f and m = %4.0f \neps = %16.8f\nNumber of Success= %4.0f\nNumber of failures=%4.0f\nAverage number of iterations = %16.4f\nAverage cpu time = %16.4f\n\n',n,m,eps(t),nbOf2_newton,nbOf1_newton,avIter_newton,av_cpu_time_newton);
        
        fprintf('FISTA :\nfor n = %4.0f and m = %4.0f \neps = %16.8f\nNumber of Success= %4.0f\nNumber of failures=%4.0f\nAverage number of iterations = %16.4f\nAverage cpu time = %16.4f\n\n',n,m,eps(t),nbOf2_fista,nbOf1_fista,avIter_fista,av_cpu_time_fista);
        
    end
end
