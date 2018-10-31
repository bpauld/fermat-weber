%   Program which compares the CPU-times used by the Weiszfeld 
%   method and by the Newton method, excluding the initialization
%   phase and based on random examples of different dimensions.

%   Some initial settings

cpu_Weiszfeld = 0;
cpu_Newton = 0;

iter_Weiszfeld = 0;
iter_Newton = 0;

fail_Weiszfeld = 0;
fail_Newton = 0;

iter_Weiszfeld = 0;
iter_Newton = 0;

number_anchor_solution = 0;

eps = 1e-5;     % parameter for termination criterion
kmax = 1000;    % maximum number of iterations for both methods
beta = 0.5;     % stepsize parameter in Newton method
sigma = 1e-4;   % stepsize parameter in Newton method
mGLL = 10;      % nonmonotone line search parameter
eta = 1e-8;     % switching criterion in Newton method
pow = 2.1;      % switching criterion in Newton method

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Outer loop for number of test runs
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

numbertestruns = 100;
counter = 0 ;

for oo = 1 : numbertestruns

oo   % output of outer iteration number

%   Generate random test problem (dimension can be modified here).
%   omega = vector of weights
%   A = matrix with the anchor points in its rows

n = 2;
m = 10;

omega = 100*rand(m,1);
A = 100*rand(m,n);





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Phase I: Initialization
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Compute minimal function value of anchor points

fp = 0;
for i = 1:m
   fp = fp + omega(i)*norm(A(1,:)-A(i,:));
end

p = 1;
fmin = fp;

for i = 2:m
   fp = 0;
   for j = 1:m
      fp = fp + omega(j)*norm(A(i,:)-A(j,:));
   end
   if fp < fmin
      p = i;
      fmin = fp;
   end
end


%   Check optimality of minimal anchor point

Rp = zeros(n,1);
for i = 1:m
   if i ~= p
      adiff = A(p,:) - A(i,:);   
      Rp = Rp + omega(i)*adiff'/norm(adiff);
   end
end
Rp = Rp';
normRp = norm(Rp,2);

stop = 1;
if normRp <= omega(p)+1e-8
   disp('Minimal anchor point is solution.')
   p;
   A(p,:)';
   stop = 0;   % to avoid using further iterations;
   number_anchor_solution = number_anchor_solution + 1;
end

%   Compute descent direction

dp = - Rp'/normRp;

%   Compute stepsize with guaranteed reduction of function value

Lp = 0;
for i = 1:m
   if i ~= p
      Lp = Lp + omega(i)/norm(A(p,:)-A(i,:),2);
   end
end

tp = (normRp - omega(p))/Lp;


%   Initial point

x0 = A(p,:)' + tp*dp;

fx0 = 0;
for i = 1:m
   fx0 = fx0 + omega(i)*norm(x0-A(i,:)');
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Phase II: Weiszfeld iteration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = x0;

%   Start to measure CPU time

tic

%   Initialization of counter

k = 0;

%   Evaluation of fixed point operator

fac = 0;
for i = 1:m
   diffx = x - A(i,:)';
   fac = fac + omega(i)/norm(diffx);
end

Tx = zeros(n,1);
gx = zeros(n,1);
for i = 1:m
   diffx = x - A(i,:)';
   normdiffx = norm(diffx,2);
   Tx = Tx + (omega(i)*A(i,:)')/normdiffx;
   gx = gx + omega(i)*diffx/normdiffx;
end
Tx = Tx/fac;

normTx = norm(x-Tx);
normgx = norm(gx,inf);

%  Initial output

disp('Initial output:')
disp('=========================================')
disp('  k      || x - T(x) ||   || grad(x) ||')
disp('=========================================')
fprintf('%4.0f %16.8f %16.8f\n',k,normTx,normgx)

while (normgx > eps) && (k < kmax) && (stop == 1)
   
   %   Weiszfeld iteration

   x = Tx;

   %   Evaluation of fixed point operator

   fac = 0;
   for i = 1:m
      diffx = x - A(i,:)';
      fac = fac + omega(i)/norm(diffx);
   end

   Tx = zeros(n,1);
   gx = zeros(n,1);
   for i = 1:m
      diffx = x - A(i,:)';
      normdiffx = norm(diffx,2);
      Tx = Tx + (omega(i)*A(i,:)')/normdiffx;
      gx = gx + omega(i)*diffx/normdiffx;
   end
   Tx = Tx/fac;

   normTx = norm(x-Tx);
   normgx = norm(gx,inf);
   k = k + 1;

   %   Output at each iteration

   fprintf('%4.0f %16.8f %16.8f\n',k,normTx,normgx);
   
end

t_Phase2_Weiszfeld = toc;
cpu_Weiszfeld = cpu_Weiszfeld + t_Phase2_Weiszfeld;
iter_Weiszfeld = iter_Weiszfeld + k;

if k>=kmax
   disp('Weiszfeld method reached maximum iteration.')
   fail_Weiszfeld = fail_Weiszfeld + 1;
end


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
   fx = fx + omega(i)*norm(x - A(i,:)');
end


%   Settings for nonmonotone line search

fvec = fx*ones(mGLL,1);
fmax = fx;
mk = 1;

%   Gradient evaluation

gx = zeros(n,1);
for i = 1:m
   diffx = x - A(i,:)';
   gx = gx + omega(i)*diffx/norm(diffx,2);
end
normgx = norm(gx,inf);



disp('Newton iteration:')
disp('============================================================')
disp('   k         f(x^k)        || Df(x^k) ||        t_k    step')
disp('============================================================')
fprintf('%4.0f %16.8f %16.8f\n',k,fx,normgx)

while (normgx > eps) && (k < kmax) && (stop == 1)
   
   %   Hessian evaluation

   Hx = zeros(n,n);
   for i = 1:m
      diffx = x - A(i,:)';
      normdiffx = norm(diffx,2);
      Hx = Hx + omega(i)*(normdiffx^2*eye(n,n)-diffx*diffx')/(normdiffx^3);
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

   %   Compute Armijo line search

   t = 1;
   xnew = x + t*d;
   fxnew = 0;
   for i = 1:m
      fxnew = fxnew + omega(i)*norm(xnew - A(i,:)');
   end
   sigmagxd = sigma*gx'*d;
   while fxnew > fmax + t*sigmagxd
      t = beta*t;
      xnew = x + t*d;
      fxnew = 0;
      for i = 1:m
         fxnew = fxnew + omega(i)*norm(xnew - A(i,:)');
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
      diffx = x - A(i,:)';
      gx = gx + omega(i)*diffx/norm(diffx,2);
   end
   normgx = norm(gx,inf);

   %   Output at each iteration

   fprintf('%4.0f %16.8f %16.8f %14.6f %4.0f\n',k,fx,normgx,t,step);
   
end

t_Phase2_Newton = toc;
cpu_Newton = cpu_Newton + t_Phase2_Newton;
iter_Newton = iter_Newton + k;

if k>=kmax
   disp('Newton method reached maximum iteration.')
   fail_Newton = fail_Newton + 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   End of outer loop
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end
x

%   Computation of average CPU-time

disp('Average CPU-time of Weiszfeld and Newton method:')
av_cpu_Weiszfeld = cpu_Weiszfeld/numbertestruns
av_cpu_Newton = cpu_Newton/numbertestruns

disp('Average number of iterations for Weiszfeld and Newton method:')
av_iter_Weiszfeld = iter_Weiszfeld/numbertestruns
av_iter_Newton = iter_Newton/numbertestruns

disp('Number of examples an anchor point was a solution:')
number_anchor_solution

disp('Number of failures of Weiszfeld and Newton iterations:')
fail_Weiszfeld
fail_Newton

counter
