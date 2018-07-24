
function [xopt, gapvalue, k] = affine(A,b,c,Mval,rhoval) 
% matrix A dimension(n,m)
% column vector b
% column vector c
% Mval: value to multiply max(abs(c)) for penalty c'x (big M-method)
% rhoval: value of step length
% Solve affine scale algorithm:  min c'x 
%                               s.t. Ax = b
%                                    x >= 0
% TEST PROBLEM
% A = [2 1 -1 0; 3 4 0 1]
% b = [2 12]
% c = [3 1 0 0]
% affine(A,b',c')

% if you use .mat LPnetlib files with bounds x.lo >= 0, x.hi <= Inf use:
% A = full(Problem.A) make double matrix
% b = Problem.b
% c = Problem.aux.c
% affine(A,b,c)

% For compare value with linprog function
% x = linprog(c,[],[],A,b,Problem.aux.lo,Problem.aux.hi);
% If you want to save more information
% [x,fval,exitflag,output]  = linprog(Problem.aux.c,[],[],full(Problem.A),Problem.b,Problem.aux.lo,Problem.aux.hi);

prompt = 'Do you want graph display? Y/N [Y]: ';
graph = input(prompt,'s');
if isempty(graph)
    graph = 'Y';
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

if(rank(A) ~= m)
     error('Error: Matrix A incorrect');
     
end
 
 optgap = 1e-5;
 myeps = 1e-12;
 M = Mval*max(abs(c));
 x = ones(n,1);
 A = [A b-A*x]; 
 % Extend A with column r
 % r = b - Aen (en is ones(n,1)
 x = [x;1]; c = [c;M];
 
 D = diag(x)^2;
 % consider identity matrix
 %D = eye(size(x,1),size(x,1));
 
 y = (A*D*A')\(A*D*c);
 k=0;
 dualgap = abs(c'*x-b'*y)/(0.1+abs(c'*x));
 gapvalue = dualgap;
 
 % bucle
 while( dualgap > optgap)
     
 z = c-A'*y; % cost - value
 dx = -D*z; % increment 
 
 if all(dx>myeps)
     x= Inf(n,1); 
     fprintf('Unbounded problem\n');
    break;
 end
     
alphaa = rhoval*min(-x(dx<myeps) ./ dx(dx<myeps));

% next x(k+1
k=k+1;
x = x+alphaa*dx;
D = diag(x)^2;
y = (A*D*A')\(A*D*c);

% compute gap
dualgap = abs(c'*x-b'*y)/(0.1+abs(c'*x));
% save gap
gapvalue = [gapvalue;dualgap];
 end 

% end means last position, so is x(n+1)
if (x(end) > 1e-4)
    fprintf('Infactible problem because need artificial variable\n');
    xopt = NaN(n,1);
    
else 
    
    fprintf('\n Value of artificial variable: %.4d\n\n',x(end))
    fprintf('Final gap value: %.2d\n\n',dualgap)
    xopt=x(1:n);
    disp('Value of vars. X')
    disp(xopt)
    
    if(graph =='Y')
    plot(0:k,gapvalue)
    title('Perform of algorithm')
    xlabel('nº iter')
    ylabel('value of gap')    
    end
end 
end

 
 
 
 
     
       
       