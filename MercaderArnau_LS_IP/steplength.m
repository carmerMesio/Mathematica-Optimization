function [alpha, alphax, alphas] = steplength(x, s, Dirx, Dirs, rho)
% 
% We need x + alphax*Dirx>0
%         s + alphas*Dirs>0

% rho indicates the maximum fraction of
% step to the boundary (value close to 1, 
% we fix in parameters rhoMin = 0.995)

% Calculate the ratio for x,s variables

  alphax = -1/min(min(Dirx./x),-1); alphax = min(1, rho * alphax);
  alphas = -1/min(min(Dirs./s),-1); alphas = min(1, rho * alphas);
  alpha = min(alphax, alphas);
  
