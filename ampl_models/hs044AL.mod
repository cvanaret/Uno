param rho_ >= 0, default 100;
param l{1..6} := 1;

var x {1..4} >= 0;
var s {1..6} >= 0;

minimize obj: 
	x[1] - x[2] - x[3] - x[1]*x[3] + x[1]*x[4] + x[2]*x[3] - x[2]*x[4]
      -l[1]*(x[1] + 2*x[2] + s[1] - 8)
      -l[2]*(4*x[1] + x[2] + s[2] - 12)
      -l[3]*(3*x[1] + 4*x[2] + s[3] - 12)
      -l[4]*(2*x[3] + x[4] + s[4] - 8)
      -l[5]*(x[3] + 2*x[4] + s[5] - 8)
      -l[6]*(x[3] + x[4] + s[6] - 5)
+ rho_/2*( (x[1] + 2*x[2] + s[1] - 8)^2 + (4*x[1] + x[2] + s[2] - 12)^2 + (3*x[1] + 4*x[2] + s[3] - 12)^2
           + (2*x[3] + x[4] + s[4] - 8)^2 + (x[3] + 2*x[4] + s[5] - 8)^2 + (x[3] + x[4] + s[6] - 5)^2 )
;

# subject to constr1: x[1] + 2*x[2] + s[1] = 8;
# subject to constr2: 4*x[1] + x[2] + s[2] = 12;
# subject to constr3: 3*x[1] + 4*x[2] + s[3] = 12;
# subject to constr4: 2*x[3] + x[4] + s[4] = 8;
# subject to constr5: x[3] + 2*x[4] + s[5] = 8;
# subject to constr6: x[3] + x[4] + s[6] = 5;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;
let x[4] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 3;
#let x[3] := 0;
#let x[4] := 4;

display obj;

solve;

display x;

display obj;

display obj + 15;
