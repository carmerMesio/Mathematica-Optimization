###Mod ampl bender
set N;
set O within N;
set A within N cross N; #arcos fijos
set Ahat within N cross N; #arcos añadibles
set AA:=A union Ahat; #todos los arcos

param xc {N};
param yc {N};
param t {i in N,l in O};
param c {(i,j) in AA, l in O} :=expresion; #arcos es decir de arco i a arco j.
param f {(i,j) in Ahat} := expresion; #costos inversio per als nous arcs que voldriem afegir (compra + instalacio)
param yb {(i,j) in Ahat}; #arcos de las nuevas posibles rutas a agregar.
param rho>0; #fluxe de entrada, va dir valor arbitrariament gran (major a la suma dels origens)
param Niter;
param nCUT;
param restric {(i,j) in Ahat, l in O,k in 1..nCUT};
param ybk {(i,j) in Ahat,k in 1..nCUT};

node I {i in N, l in O}: net_out=t[i,l]; # si positivo ==> inyección, 
                                         # si negativo extracción
arc xl {(i,j) in AA, l in O}>=0: from I [i,l], to I [j,l]; ##capacidad de arco i,j. Valor que pasa d'un punt a l'altre del arc.

## node I i xl es la restriccio (1) Bxl = tl xl>=0.
## net out fa referencia a on pot anar el flux o la capacitat maxima del node. es fa servir net_out perque ampl entengui que tot el que surt del nodo origen O ha d'acabar sen 0
## es a dir si mirem a t el que surt de O que es 100 es compensa amb tots els altres elements de la fila O.

var y{(i,j) in Ahat} binary; ##para ver si se activa o no. Del conjunt d'arcs que volem afegir si el fem servir o no.

#Problema original
minimize z: #para cada par de nodos en Ahat calculamos su f (costes de inversion) y los costes para cada arco.
   sum {(i,j) in Ahat} f[i,j]*y[i,j]+
   sum{l in O} (sum {(i,j) in AA} c[i,j,l]*xl[i,j,l]);

subject to caps1 {(i,j) in Ahat, l in O}:
    xl[i,j,l]<=rho*y[i,j];

#Subproblema 
minimize zd: sum{l in O} (sum {(i,j) in AA} c[i,j,l]*xl[i,j,l]);

subject to caps {(i,j) in Ahat, l in O}:
    xl[i,j,l]<=rho*yb[i,j];

#Master problem (yb)
#param u {i in N, l in O,k in 1..nCUT}<=0;
#var zmp;

#minimize ZMP:zmp;

#
#subject to Bcut {k in 1..nCUT}:
#zmp>=(sum {(i,j) in Ahat} f[i,j]*y[i,j])+
#   (sum{i in N, l in O} t[i,l]*u[i,l,k])+
 #  rho*(sum{(i,j) in Ahat, l in O} restric[i,j,l,k]*(1-ybk[i,j,k])*y[i,j]);