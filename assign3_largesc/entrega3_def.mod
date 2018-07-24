set nusos;										
set centr within nusos;							
set links within ( nusos cross nusos );			
set origens within centr;						
set destins within centr;						
set odpair within ( origens cross destins );	

# Sets of destinations by origin
set destperorig { i in origens } :=	setof { (i1, j1) in odpair : i = i1 } j1;

param g{odpair} > 0;							

# Cost of the arc function: t[i,j] = c[i,j] + d[i,j] * x
param d{links};
param c{links};									
param t0{links};								
							# Capacity of the arcs
param rho default 2;

## Define working sets for the comparision.
param Wx{links} default 0;
param Ws{links, 1..rho} default 0;
param W{links, 0..rho} default 0;
param flag{0..rho}>=0; #save the positions ocupped 


# Right side coefficient of the balance equation
param Tdreta { i in nusos, k in origens }:=

if i in destperorig[k] then -1.0*g[k, i]
else 
if i = k then sum {j in destperorig[k]} g[k, j]
else 0;
# Right side coefficient of the balance equation

# Balance equation
node N {i in nusos, k in origens}: net_out = Tdreta[i, k];

# Variables
arc v_k { (i, j) in links, k in origens } >= 0, from N[i, k], to N[j, k] ;	
var v { (i, j) in links };													
var vvv{ (i, j) in links };
var alp{l in 0..rho} >= 0;

####### SUBPROBLEM (LINEARIZATION)
minimize Vg: sum { (i, j) in links } t0[i, j]*v[i, j];

# Flow definition constraint
subject to flux_total { (i, j) in links }:
v[i, j] = sum { k in origens } v_k[i, j, k];


####### MASTER (QUADRATIC FUNCTION)

minimize Vnl: sum {(i, j) in links }c[i,j]*vvv[i,j] + sum {(i, j) in links }0.5*d[i,j]*vvv[i,j]^2; 

subject to suma: sum{l in 0..rho}alp[l] = 1;

subject to wset {(i,j) in links}: 
vvv[i,j]= sum{l in 0..rho}W[i,j,l]*alp[l];

