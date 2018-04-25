var x {1..4} >= -10, <= 10;

minimize obj: 
  100*(x[2]-x[1]^2)^2 + (1-x[1])^2 + 90*(x[4]-x[3]^2)^2 + (1-x[3])^2
  + 10.1*( (x[2]-1)^2 + (x[4]-1)^2 ) + 19.8*(x[2]-1)*(x[4]-1)
  ;

#subject to constr1: x[1] + 2*x[2] + 2*x[3] <= 72;
#subject to constr2: x[1] + 2*x[2] + 2*x[3] >= 0;

let x[1] := -3;
let x[2] := -1;
let x[3] := -3;
let x[4] := -1;

#printf "optimal solution as starting point \n";
#let x[1] := 1;
#let x[2] := 1;
#let x[3] := 1;
#let x[4] := 1;

display obj;

solve;

display x;

display obj;

display obj - 0;
