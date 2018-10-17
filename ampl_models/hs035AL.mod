param rho_ >= 0, default 10000;
param l := 0.25;

var x {1..3} >= 0;
var s >= -3;

minimize obj: 
  9 - 8*x[1] - 6*x[2] - 4*x[3] + 2*x[1]^2 + 2*x[2]^2 + x[3]^2 
  + 2*x[1]*x[2] + 2*x[1]*x[3]
  - l*(x[1] + x[2] + 2*x[3] + s)
  + rho_/2*(x[1] + x[2] + 2*x[3] + s)^2
;

##subject to constr1: x[1] + x[2] + 2*x[3] + s = 3;

let x[1] := 0.5;
let x[2] := 0.5;
let x[3] := 0.5;

#printf "optimal solution as starting point \n";
#let x[1] := 4/3;
#let x[2] := 7/9;
#let x[3] := 4/9;

display obj;
expand;

solve;

display x,s,l;

display (x[1] + x[2] + 2*x[3] + s)^2;

display obj;

display obj - 1/9;
