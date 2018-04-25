param a {1..2, 1..4};
param b {1..2};

var x {1..4} >= 0.001;

minimize obj: 
  1 + sum {j in 1..4} x[j]
  ;

subject to constr {i in 1..2}: sum {j in 1..4} a[i,j]/x[j] <= b[i];
subject to ub {j in 1..4}: x[j] <= (5-j)*1e5;

let x[1] := 1;
let x[2] := 1;
let x[3] := 1;
let x[4] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 193.4071;
#let x[2] := 179.5475;
#let x[3] := 185.0186;
#let x[4] := 168.7062;

data;

param a: 1 2 3 4 :=
    1  4    2.25 1    0.25
    2  0.16 0.36 0.64 0.64
    ;

param b :=
    1 0.0401
    2 0.010085
    ;

display obj;

solve;

display x;

display obj;

display obj - 727.67937;
