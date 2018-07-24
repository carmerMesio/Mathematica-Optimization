function [x,y,f,hist] = path_corrected(A,b,c,option)
% Algorithm Primal/Dual IP
% Description of problem:
% PRIMAL: min c'x s.t. Ax=b, x>=0, 
% DUAL:   max b'y  s.t. A'y+s=c,  s>=0%
%
% input: A (mxn) matrix (in format sparse to apply cholesky)
%        b rhs vector
%        c cost vector
%        
%        option  if is >= 1 perform heuristic Mehrotra, else perform 
%        primal dual with constant sigma = 0.1
%        by deffect: algorithm perform Mehrotra
%
% output: x primal vars solution
%         y dual vars solution
%         f optimal objective value c'x = b'y (because of LP)
%         
% internal parameters: 
%         itmax max iterations
%         tol convergence tolerance
%         maxDiag element for the X^{-1}S matrix to avoid errors
%         rhoMin minimum value of the step length

% check validity of input arguments
if (nargin < 3)
  error('Put A,b,c please...');
end

% Peform Mehrotra by deffect
if isempty(option)
    option = 1;
end


% put help issparse in Matlab to reference
if issparse(A) < 1
  A = sparse(A);
end

[m,n] = size(A);
[mc,nc] = size(c);
[mb,nb] = size(b);

% ~= means !=
if (or(nc ~=1,nb ~=1))
    error('Error: b or/and c are not column vectors (n,m=1)');
end

if (or(mc ~= n,mb ~=m))
    error('Error: dimension incorrect');
end

prompt = 'Do you want graph display? Y/N [Y]: ';
graph = input(prompt,'s');
	if isempty(graph)
		graph = 'Y';
	end


% empty plots for output display
hist_k = [];
hist_f = [];
hist_mu = [];
hist_nres = [];

% Internal initial parameters
itmax = 200; 
tol = 1.e-7;
maxDiag = 5.e+15;
rhoMin = 0.995;

% start the clock
t0=cputime;


% set initial point, max(abs(c))
bigM = 100*max(abs(c));
x = bigM*ones(n,1); s = x; y = zeros(m,1);


% find row/column ordering for a sparse Cholesky factorization of ADA'
ordering = symamd(A*A');
bc = 1+max([norm(b), norm(c)]);

for iter=1:itmax
  
% compute residuals
  Rd = A'*y+s-c;  % dual errors
  Rp = A*x-b; % primal errors
  Rc = x.*s; % complementary errors
  mu = mean(Rc); % update mu (proximity to optimallity)
   
 
% check relative decrease in residual, for purposes of convergence test
  relResidual = norm([Rd;Rp;Rc])/bc;
  fprintf(1,'iter %2i: mu = %9.2e, resid = %9.2e\n', iter, full(mu), ...
	  full(relResidual));

% test for convergence and break
  if(relResidual <= tol & mu <= tol) break; end;

 if option >= 1 
 
  % Mehrotra heuristic approach
  
  % set up the scaling matrix and form the coef matrix for normal equations
  d = min(maxDiag, x./s);
  B = A*sparse(1:n,1:n,d)*A';
  
  % Cholesky with permuting
  R = cholinc(B(ordering,ordering),'inf');
  
  % set up the right-hand side
  t1 = x.*Rd-Rc;
  t2 = -(Rp+A*(t1./s));
  
  
  % solve the normal equations system for dy and recover the other 
  %  step components dx and ds
  dy = zeros(m,1);
  dy(ordering) = R\(R'\t2(ordering));
  dx = (x.*(A'*dy)+t1)./s;
  ds = -(s.*dx+Rc)./x;
  
  
  
% set the parameter eta defining fraction of max step to boundary
% rho = max(rhoMin,1-mu);
[alpha, alphax, alphas] = steplength(x, s, dx, ds, 1);
  
% Central parameter SIGMA with heuristic
%mu_a = mean(((x + alphax*dx).*(s + alphas*ds)));
mu_a = (x+alphax*dx)'*(s+alphas*ds)/n;
sigma = (mu_a/mu)^3;

% Change complementary restriction
Rc = Rc + dx.*ds - sigma*mu*ones(n,1);

% CORRECTOR STEP NOW
% repeat normal equations with same Cholesky matrix R

t1 = x.*Rd-Rc;
t2 = -Rp -A*(t1./s);

%t2 = -Rp - A*(d.*Rd-Rc./s);
dy1 = zeros(m,1);
dy1(ordering) = R\(R'\t2(ordering));
dx1 = (x.*(A'*dy1)+t1)./s;
ds1 = -(s.*dx1+Rc)./x;

% Combine two directions
dy1 = dy1 +  dy;
dx1 = dx1 +  dx;
ds1 = ds1 +  ds;

rho = max(rhoMin,1-mu);
[alpha, alphax, alphas] = steplength(x, s, dx1, ds1, rho);
  
% take the step
x = x + alphax * dx1;
s = s + alphas * ds1;
y = y + alphas * dy1;

f = c'*x;
hist_f = [hist_f;f];
hist_k = [hist_k;iter];
hist_mu = [hist_mu;mu];
hist_nres = [hist_nres;relResidual];

 else
   
 % Algorithm normal with effect of KKT-mu with constant sigma 0.1
 sigma = 0.1;
 Rc = Rc - sigma*mu;
  
  % set up the scaling matrix and form the coef matrix for normal equations
  d = min(maxDiag, x./s);
  B = A*sparse(1:n,1:n,d)*A';
  
  % Cholesky with permuting
  R = cholinc(B(ordering,ordering),'inf');
  
  % set up the right-hand side
  t1 = x.*Rd-Rc;
  t2 = -(Rp+A*(t1./s));
  
  
% solve the normal equations system for dy and recover the other 
%  step components dx and ds
dy = zeros(m,1);
dy(ordering) = R\(R'\t2(ordering));
dx = (x.*(A'*dy)+t1)./s;
ds = -(s.*dx+Rc)./x;
  
  
% resize rho if are nearly to the opt solution
rho = max(rhoMin,1-mu);
[alpha, alphax, alphas] = steplength(x, s, dx, ds, rho);
  
% take the step
x = x + alphax * dx;
s = s + alphas * ds;
y = y + alphas * dy;

f = c'*x;

hist_f = [hist_f;f];
hist_k = [hist_k;iter];
hist_mu = [hist_mu;mu];
hist_nres = [hist_nres;relResidual];


     end
end


% >200 iterations
if(or(relResidual > tol,mu > tol))
    fprintf('No solution found\n');
    % save list of results (output param)
    hist = [hist_k,hist_f,hist_mu,hist_nres];


else
    
% calculate value of f.obj (remember c'x = b'y)
f = c'*x;
x=full(x); s=full(s); y=full(y);

% save list of results (output param)
hist_f = [hist_f;f];
hist_k = [hist_k;iter];
hist_mu = [hist_mu;mu];
hist_nres = [hist_nres;relResidual];
hist = [hist_k,hist_f,hist_mu,hist_nres];

fprintf('Solution found!\t CPU time required = %g\n', cputime-t0);
fprintf('Value of f.obj %9.4e\n', f);


end

if(graph =='Y')
% save all values
k = hist(:,1);
ff = hist(:,2);
mu = hist(:,3);
res = hist(:,4);


subplot(1,3,1)
set(gca,'FontSize',13)
% semilogy plots data with logarithmic scale for the y-axis.
plot(k(1:end),ff(1:end),'-rs','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',2)
xlabel('# Iteration')
ylabel('Value f.obj')
subplot(1,3,2)
set(gca,'FontSize',13)
% % semilogy plots data with logarithmic scale for the y-axis.
semilogy(k(1:end),mu(1:end),'-rs','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',2)
xlabel('# Iteration')
ylabel('log(\mu)')
title( {'Performance IP algorithm';sprintf('Value of f.obj = %g',f)})
subplot(1,3,3)
set(gca,'FontSize',13)
semilogy(k(1:end),res(1:end),'-rs','LineWidth',2,'MarkerEdgeColor','k','MarkerSize',2)
xlabel('# Iteration')
ylabel('log(rel residual)')


end



return;  
