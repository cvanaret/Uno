param l{1..6};
param u{1..6};

param rho := 0.2;
param mu{1..6};
param sigma{1..6};

var y {j in 1..6} >= (l[j] - mu[j])/sigma[j], <= (u[j]-mu[j])/sigma[j];

minimize obj: 
      ( y[1]^2 + 2*rho*y[1]*y[2] + y[2]^2 ) / (1-rho^2)
      + sum {j in 3..6} y[j]^2
  ;

subject to constr1: 
	y[1]/sigma[2] + 4000*y[2]/sigma[1] = 2000/sigma[1] + 0.2/sigma[2];

data;

param l :=
  1    0
  2  -10
  3    0
  4    0
  5   -1
  6    0
  ;

param u :=
  1  2.0e+4
  2  10
  3  1.0e+7
  4  20
  5  1
  6  2.0e+8
  ;

param mu :=
  1  10000
  2      1
  3  2e+6
  4     10
  5      0.001
  6  1e+8
  ;

param sigma :=
  1  8000
  2      1
  3  7e+6
  4     50
  5      0.05
  6  5e+8
  ;

let y[1] := 6.0e+3;
let y[2] := 1.5;
let y[3] := 4.0e+6;
let y[4] := 2;
let y[5] := 3.0e-3;
let y[6] := 5.0e+7;

#printf "optimal solution as starting point \n";
#let y[1] := 91600/7;
#let y[2] := 79/70;
#let y[3] := 2.0e+6;
#let y[4] := 10;
#let y[5] := 1.0e-3;
#let y[6] := 1.0e+8;

let {j in 1..6} y[j] := (y[j] - mu[j])/sigma[j];

display -exp(-obj/2);

option loqo_options $loqo_options" convex";

solve;

display -exp(-obj/2);

display -exp(-obj/2) + exp(-27/280);
