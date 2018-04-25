param u {1..6};

var x {j in 1..6} >= 0, <= u[j];

minimize obj: 
  4.3*x[1] + 31.8*x[2] + 63.3*x[3] + 15.8*x[4] + 68.5*x[5] + 4.7*x[6]
  ;

subject to constr1: 
    17.1*x[1] + 38.2*x[2] + 204.2*x[3] + 212.3*x[4] + 623.4*x[5] + 1495.5*x[6]
    - 169*x[1]*x[3] - 3580*x[3]*x[5] - 3810*x[4]*x[5] - 18500*x[4]*x[6]
    - 24300*x[5]*x[6] >= 4.97;

subject to constr2: 
    17.9*x[1] + 36.8*x[2] + 113.9*x[3] + 169.7*x[4] + 337.8*x[5] + 1385.2*x[6]
    - 139*x[1]*x[3] - 2450*x[4]*x[5] - 16600*x[4]*x[6] - 17200*x[5]*x[6]
    >= -1.88;

subject to constr3: 
    -273*x[2] - 70*x[4] - 819*x[5] 
    + 26000*x[4]*x[5] >= -29.08;

subject to constr4: 
    159.9*x[1] - 311*x[2] + 587*x[4] + 391*x[5] + 2198*x[6]
    - 14000*x[1]*x[6] >= -78.02;

data;

param u :=
  1   0.31
  2   0.046
  3   0.068
  4   0.042
  5   0.028
  6   0.0134
  ;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;
let x[4] := 0;
let x[5] := 0;
let x[6] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 0;
#let x[3] := 0;
#let x[4] := 0;
#let x[5] := 0;
#let x[6] := 0.0033233033;

solve;

display x;

display obj;

display obj - 0.015619514;
