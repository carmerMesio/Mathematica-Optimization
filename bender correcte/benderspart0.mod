set N;
set O within N;
set A within N cross N; #arcos fijos
set Ahat within N cross N; #arcos añadibles
set AA:=A union Ahat; #todos los arcos

param xc {N};
param yc {N};
param t {i in N,l in O};
param c {(i,j) in AA, l in O} := 100+5*sqrt((xc[i]-xc[j])**2+(yc[i]-yc[j])**2);
param f {(i,j) in Ahat} := 22*(2*abs(xc[i]-xc[j])+1.5*abs(yc[i]-yc[j]));
param yb {(i,j) in Ahat};
param rho>0;

node I {i in N, l in O}: net_out=t[i,l]; # si positivo ==> inyección, 
                                         # si negativo extracción
arc xl {(i,j) in AA, l in O}>=0: from I [i,l], to I [j,l];

var y{(i,j) in Ahat} binary;


#Problema original
minimize z: 
   sum {(i,j) in Ahat} f[i,j]*y[i,j]+
   sum{l in O} (sum {(i,j) in AA} c[i,j,l]*xl[i,j,l]); #cost associat a portar 10 unitats de fluxe de i a j.

subject to caps1 {(i,j) in Ahat, l in O}:
    xl[i,j,l]<=rho*y[i,j];

	######### Parte 1 hasta aqui###########