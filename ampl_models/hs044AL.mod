param rho_ >= 0, default 1000;
param l{1..6}, default 0;
let l[3] := 1.25;
let l[5] := 1.5;
param sl{1..6} default -8;
let sl[2] := -12;
let sl[3] := -12;
let sl[6] := -5;

var x {1..4} >= 0;
var s {i in 1..6} >= sl[i];

minimize obj: 
	x[1] - x[2] - x[3] - x[1]*x[3] + x[1]*x[4] + x[2]*x[3] - x[2]*x[4]
      -l[1]*(x[1] + 2*x[2] + s[1])
      -l[2]*(4*x[1] + x[2] + s[2])
      -l[3]*(3*x[1] + 4*x[2] + s[3])
      -l[4]*(2*x[3] + x[4] + s[4])
      -l[5]*(x[3] + 2*x[4] + s[5])
      -l[6]*(x[3] + x[4] + s[6])
+ rho_*( (x[1] + 2*x[2] + s[1])^2 + (4*x[1] + x[2] + s[2])^2 + (3*x[1] + 4*x[2] + s[3])^2
           + (2*x[3] + x[4] + s[4])^2 + (x[3] + 2*x[4] + s[5])^2 + (x[3] + x[4] + s[6])^2 )
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
display (x[1] + 2*x[2] + s[1])^2 + (4*x[1] + x[2] + s[2])^2 + (3*x[1] + 4*x[2] + s[3])^2
           + (2*x[3] + x[4] + s[4])^2 + (x[3] + 2*x[4] + s[5])^2 + (x[3] + x[4] + s[6])^2 ;

display obj + 15;
