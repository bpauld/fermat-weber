function plotBFGSvsSub(A,omega,nbOfValues)
x_axis = 1:nbOfValues ;
[n,m] = size(A) ;


x0 = -ones(n,1) ;
c1 = 0.25 ;
c2 = 0.5 ;
eps1 = 1e-8 ;
eps2 = 1e-6 ;

[x,k,stop] = bfgs_inf(A,omega,x0,c1,c2,eps1,eps2) ;

f_bfgs = zeros(nbOfValues,1) ;
first_f_bfgs = bfgs_inf_values(A,omega,x0,c1,c2,k) ;
for i=1:k
    f_bfgs(i,1) = first_f_bfgs(i,1) ;
end
for i=k+1:nbOfValues
    f_bfgs(i,1) = first_f_bfgs(k,1) ;
end

f_sub = subgradient_values(A,omega,inf,nbOfValues) ;

figure
%plot(x_axis,f_bfgs,'--',x_axis,f_sub) ;
semilogx(x_axis,f_bfgs,'--',x_axis,f_sub) ;
xlabel("Iteration")
ylabel('Function value')
title("Convergence for l-infinity norm")
legend('BFGS','Subgradient Method')

%print -deps bfgsVSsubLOG %to save figure

