var x {1..3} >= 0;

minimize obj: 
  9 - 8*x[1] - 6*x[2] - 4*x[3] + 2*x[1]^2 + 2*x[2]^2 + x[3]^2 
  + 2*x[1]*x[2] + 2*x[1]*x[3]
  ;

subject to constr1: x[1] + x[2] + 2*x[3] <= 3;

let x[1] := 0.5;
let x[2] := 0.5;
let x[3] := 0.5;

#printf "optimal solution as starting point \n";
#let x[1] := 4/3;
#let x[2] := 7/9;
#let x[3] := 4/9;

display obj;

solve;

display x;

display obj;

display obj - 1/9;
