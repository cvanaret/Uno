param a;
param l {1..4};
param u {1..4};

var x {j in 1..4} >= l[j], <= u[j];

minimize obj: 
  3*x[1] + 1.0e-6*x[1]^3 + 2*x[2] + 2.0e-6*x[2]^3/3
  ;

subject to constr1: -a <= x[4] - x[3] <= a;
subject to constr2: x[1] = 1000*sin(-x[3] - 0.25) +
			   1000*sin(-x[4] - 0.25) +
			   894.8;
subject to constr3: x[2] = 1000*sin(x[3] - 0.25) +
			   1000*sin(x[3]-x[4] - 0.25) +
			   894.8;
subject to constr4: 1000*sin(x[4] - 0.25) + 1000*sin(x[4] - x[3] - 0.25) +
			   1294.8 = 0;

data;

param a := 0.55;
param l :=
  1  0
  2  0
  ;
param u :=
  1  1200
  2  1200
  ;

let l[3] := -a;
let l[4] := -a;

let u[3] :=  a;
let u[4] :=  a;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;
let x[4] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 679.9453;
#let x[2] := 1026.067;
#let x[3] :=  0.1188764;
#let x[4] := -0.3962336;

data;

solve;

display x;

display obj;

display obj - 5126.4981;
