# Model modificat per cutting planes
#per el cutting plane necessitem una solucio factible (usem variables artificials)
set N;
set artN;    #conjunt unitari
set O within N;
set A within N cross N; #arcos fijos
#set Aart1{O}; 
param M;
param y{A} default M;

param xc {N};
param yc {N};
param t {i in N,l in O};
param c {(i,j) in A} :=150+(((xc[i]-xc[j])^2)+((yc[i]-yc[j])^2)); ##idem a benders decomposition
param rho>0;
param nCUT;
param mu{(i,j) in A} >=0, default 0; # els multiplicadors
var xx{A};
##nodes original in the network
node I {i in N, l in O}: net_out=t[i,l]; #si positivo ==> inyección, si negativo extracción
node art{o in artN, l in O}: net_out=0; #links artificials no hi ha extraccio per costos elevats

#dos tipus de artificial nodes de origen a node aritificial (1) i de node artificial a destí.
arc xl {(i,j) in A, l in O}>=0: from I [i,l], to I [j,l];
arc art1{o in artN, l in O}>=0: from I [l,l], to art [o,l];
arc art2{o in artN, i in N, l in O} >=0: from art[o,l], to I[i,l];

subject to total_flow {(i,j) in A}: xx[i,j] = sum{l in O} xl[i,j,l]; ##capacitats que passen per nodes determinats.

minimize w: (sum {(i,j) in A} c[i,j]*xx[i,j])+
                sum{(i,j) in A} mu[i,j]*(xx[i,j] - y[i,j])+
				#cost terms for the artificial nodes.
                rho*(sum{o in artN, l in O} art1[o,l])+          
                rho*(sum{o in artN, i in N, l in O} art2[o,i,l]);
subject to caps {(i,j) in A}: xx[i,j]<=y[i,j];

var z;
#mu0 optimization variables 
var mu0{A}>=0;
param YY{1..nCUT} default M;
minimize Z: z;
param xxX{A,{1..nCUT}} default 0;
param art1X{o in artN,l in O, k in {1..nCUT}} default 0;
param art2X{o in artN, i in N, l in O, k in {1..nCUT}} default 0;
##cuts for solving the cutting plane algorithm. xxX original flows and then art flows
subject to cuts{k in {1..nCUT}}: z <= (sum {(i,j) in A} c[i,j]*xxX[i,j,k])+
                sum{(i,j) in A} mu0[i,j]*(xxX[i,j,k] - y[i,j])+
                rho*(sum{o in artN, l in O} art1X[o,l,k])+
                rho*(sum{o in artN, i in N, l in O} art2X[o,i,l,k]) + M*YY[k];


