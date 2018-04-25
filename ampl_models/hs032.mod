var x {1..3} >= 0;

minimize obj: 
  (x[1] + 3*x[2] + x[3])^2 + 4*(x[1] - x[2])^2
  ;

subject to constr1: 6*x[2] + 4*x[3] - x[1]^3 >= 3;
subject to constr2: x[1] + x[2] + x[3] = 1;

let x[1] := 0.1;
let x[2] := 0.7;
let x[3] := 0.2;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 0;
#let x[3] := 1;

display obj;

solve;

display x;

display obj;

display obj - 1;
