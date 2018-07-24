
function [xopt, gapvalue,i] = affine(A,b,c) 
% clear * , clear all (erase workspace)
% joc de proves
% A = [2 1 -1 0; 3 4 0 1]
% b = [2 12]
% c = [3 1 0 0]
% [bla bla] = affine(A,b',c') ... return multiple outputs ... 
% put ; for dont display  ans values...

% plot(0:i,gapvalue)
% title('Progress of algorithm')
% xlabel('nº iter')
% ylabel('gap value')



% carregar prova mat
% atacar a struct .. Problem.A , Problem.b , etc...
% hauria de donar ..


% A = full(Problem.A) per portar sparse to matrix full ?
% després b = Problem.b
% després c = Problem.aux.c

% PER PROVAR AMB LINPROG;
% x = linprog(c,[],[],A,b,Problem.aux.lo,Problem.aux.hi);
% per veure més info
% [x,fval,exitflag,output]  = linprog(Problem.aux.c,[],[],full(Problem.A),Problem.b,Problem.aux.lo,Problem.aux.hi);

[m,n] = size(A);
[mc,nc] = size(c);
[mb,nb] = size(b);

% ~= means !=
if (or(nc ~=1,nb ~=1))
    error('Error: c,b deben ser vector columna');
end 
if (or(mc ~= n,mb ~=m))
    error('Error: dimensiones incompatibles para c,b, A');
end

if(rank(A) ~= m)
     error('Error: Matriz A rango completo');
     
end
 
 optgap = 1e-6;
 myeps = 1e-12;
 rho = 0.85 ;
 M = 100*max(abs(c));
 x = ones(n,1);
 A = [A b-A*x]; 
 % Extend A with column r
 % r = b - Aen (en is ones(n,1)
 x = [x;1];
 c = [c;M];
 % D = diag(x)^2;
 

 % consider identity matrix
 D = eye(size(x,1),size(x,1));
 
 
 % buscar resoldre y sense fer inv()
 % y =inv(A*D*A')*A*D*c;
 y = (A*D*A')\(A*D*c);
 i=0;
 dualgap = abs(c'*x-b'*y)/(0.1+abs(c'*x));
 gapvalue = dualgap;
 % bucle
 while( dualgap > optgap)
     
 z = c-A'*y; % cost - resposta trobada
 dx = -D*z; % increment 
 
 if all(dx>myeps)
     x= Inf(n,1); % alternative Inf(n,1)
     fprintf('Unbounded problem\n');
     break;
 end
     
alphaa = rho*min(-x(dx<myeps) ./ dx(dx<myeps));

% proxima x(k+1
i=i+1;
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
    fprintf('Infactible problem\n');
    xopt = NaN(n,1);
    
else 
    
    fprintf('\nValor de variable artificial: %.4d\n\n',x(end))
    fprintf('final gap value: %.2d\n\n',dualgap)
    xopt=x(1:n);
    disp('Valor variables X')
    disp(xopt)
    
    
end 
end

 
 
 
 
     
       
       