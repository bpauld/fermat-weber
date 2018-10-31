function generalizedFW(H,numbertestruns,m)
%param H : symmetric positive definite matrix defining the distance
%param numbertestruns : number of test runs
%param m: number of anchor points
%Program that implements the Newton for the generalized Fermat Weber
%Problem


H

cpu_Newton = 0;
fail_Newton = 0;
iter_Newton = 0;

number_anchor_solution = 0;

eps = 1e-8;     % parameter for termination criterion
kmax = 1000;    % maximum number of iterations for both methods
beta = 0.5;     % stepsize parameter in Newton method
sigma = 1e-4;   % stepsize parameter in Newton method
mGLL = 10;      % nonmonotone line search parameter
eta = 1e-8;     % switching criterion in Newton method
pow = 2.1;      % switching criterion in Newton method

%first let's find the square root of the matrix H, which will be very
%useful
S = sqrtm(H) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Outer loop for number of test runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




for oo = 1 : numbertestruns

%oo   % output of outer iteration number

%   Generate random test problem (dimension can be modified here).
%   omega = vector of weights
%   A = matrix with the anchor points in its rows

[n,w] = size(H); %disregard w

omega = 100*rand(m,1);
A = 100*rand(m,n) ;




% compute minimum function value of anchor points

fp = 0 ;
for i=1:m
    fp = fp + omega(i)*norm(S*(A(1,:)-A(i,:))') ;
end

p = 1;
fmin = fp;

for i = 2:m
   fp = 0;
   for j = 1:m
      fp = fp + omega(j)*norm(S*(A(j,:)-A(i,:))');
   end
   if fp < fmin
      p = i;
      fmin = fp;
   end
end



%   Check optimality of minimal anchor point

crit = zeros(n,1) ;
for i=1:m
    if i~=p
        num = S*(A(p,:) - A(i,:))' ;
        crit = crit + omega(i)*num/norm(num) ;
    end
end
normCrit = norm(crit) ;


stop = 1;
if normCrit <= omega(p)
   %disp('Minimal anchor point is solution.')
   p;
   A(p,:)';
   stop = 0;   % to avoid using further iterations;
   number_anchor_solution = number_anchor_solution + 1;
   continue
end


%   Compute descent direction
Rp = zeros(n,1);
for i = 1:m
   if i ~= p
      adiff = (A(p,:) - A(i,:))';   
      Rp = Rp + omega(i)*adiff/norm(S*(A(p,:) - A(i,:))');
   end
end
normRp = norm(S*Rp);

dp = -Rp / normRp ;

%backtracking===================================================================================================
%FOR NOW WE USE THE SAME AS PAPER
%   Compute stepsize with guaranteed reduction of function value


x = A(p,:) ;



t = 2 ;
xnew = x + t*dp' ;



fxnew = 0 ;
for i =1:m
    fxnew = fxnew + omega(i)*norm(S*(xnew-A(i,:))') ;
end


while fxnew >= fmin
    t = t/2 ;
    xnew = x +t*dp' ;
    fxnew = 0 ;
    for i =1:m
        fxnew = fxnew + omega(i)*norm(S*(xnew-A(i,:))') ;
    end
end



x0 = xnew' ;


%}



%{
Lp = 0;
for i = 1:m
   if i ~= p
      Lp = Lp + omega(i)/norm(H*(A(p,:)-A(i,:))',2);
   end
end

tp = (normRp - omega(p))/Lp;



%   Initial point

x0 = A(p,:)' + tp*dp;

%}




%====================================================================================================







%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Phase II: Newton iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x0;

%   Start to measure CPU time

tic

%   Initialization of counter

k = 0;

%   Function evaluation

fx = 0;
for i = 1:m
   fx = fx + omega(i)*norm(S*(x - A(i,:)'));
end


%   Settings for nonmonotone line search

fvec = fx*ones(mGLL,1);
fmax = fx;
mk = 1;

%   Gradient evaluation

%{
gx = zeros(n,1);
for i = 1:m
   diffx = H*(x - A(i,:)');
   gx = gx + omega(i)*diffx/norm(S*(x-A(i,:)'));
end
%}
gx = l2_gradient(H,A',omega,x) ;
normgx = norm(gx,inf);

%disp('Newton iteration:')
%disp('============================================================')
%disp('   k         f(x^k)        || Df(x^k) ||        t_k    step')
%disp('============================================================')
%fprintf('%4.0f %16.8f %16.8f\n',k,fx,normgx)

while (normgx > eps) && (k < kmax) && (stop == 1)
   
   %   Hessian evaluation
   Hx = zeros(n,n);
   for i = 1:m
      diffx = x - A(i,:)';
      normdiffx = norm(S*diffx);
      Hx = Hx + omega(i)*(normdiffx^2*H-H*(diffx*diffx')*H)/(normdiffx^3);
   end
   
   %   Compute search direction
   %   step = 0: Newton direction
   %   step = 1: gradient direction (or Weiszfeld!?)

   
   condi=rcond(Hx);
   if condi>1e-6
      step = 0;
      %L = chol(Hx, 'lower');
      %dd = L\gx;
      %d = -L'\dd;
      d = - Hx\gx;  
      if (gx'*d) > -(eta*norm(d)^pow)
         d = - gx;
         step = 1;
      end
   else
      d = - gx;
      step = 1;
   end
   
  
   d = -Hx\gx ;
   
   %   Compute Armijo line search

   t = 1;
   xnew = x + t*d;
   fxnew = 0;
   for i = 1:m
      fxnew = fxnew + omega(i)*norm(S*(xnew - A(i,:)'));
   end
   %fx = 0 ;
   %for i =1:m
    %   fx = fx + omega(i)*norm(S*(x-A(i,:)')) ;
   %end
   sigmagxd = sigma*gx'*d;
   while fxnew > fmax + t*sigmagxd
      t = beta*t;
      xnew = x + t*d;
      fxnew = 0;
      for i = 1:m
         fxnew = fxnew + omega(i)*norm(S*(xnew - A(i,:)'));
      end
   end

   
   %   Updates

   x = xnew;
   fx = fxnew;
   k = k + 1;
   
   
   %   Updates nonmonotone line search

   if step == 0
      mk = min(mk+1,mGLL);
   else
      mk = 1;
   end
   fvec(2:mGLL) = fvec(1:mGLL-1);
   fvec(1) = fx;
   fmax = max(fvec(1:mk));
  
   
   %   Gradient evaluation

   gx = zeros(n,1);
   for i = 1:m
       diffx = H*(x - A(i,:)');
       gx = gx + omega(i)*diffx/norm(S*(x-A(i,:)'));
   end
   normgx = norm(gx,2);
   
   %   Output at each iteration

   %fprintf('%4.0f %16.8f %16.8f %14.6f %4.0f\n',k,fx,normgx,t,step);  
   
end

t_Phase2_Newton = toc;
cpu_Newton = cpu_Newton + t_Phase2_Newton;
if k < kmax
    iter_Newton = iter_Newton + k;
end

if k>=kmax
   %disp('Newton method reached maximum iteration.')
   fail_Newton = fail_Newton + 1;
end

end

disp('Average CPU-time of Newton method:')
av_cpu_Newton = cpu_Newton/(numbertestruns-number_anchor_solution)

disp('Average number of iterations for Newton method:')
av_iter_Newton = iter_Newton/(numbertestruns-number_anchor_solution)

disp('Number of examples an anchor point was a solution:')
number_anchor_solution

disp('Number of failures of Newton iterations:')
fail_Newton


end