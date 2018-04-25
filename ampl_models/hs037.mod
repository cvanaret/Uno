var x {1..3} >= 0, <= 42;

minimize obj: 
  -x[1]*x[2]*x[3]
  ;

subject to constr1: x[1] + 2*x[2] + 2*x[3] <= 72;
subject to constr2: x[1] + 2*x[2] + 2*x[3] >= 0;

let x[1] := 10;
let x[2] := 10;
let x[3] := 10;

#printf "optimal solution as starting point \n";
#let x[1] := 24;
#let x[2] := 12;
#let x[3] := 12;

display obj;

solve;

display x;

display obj;

display obj + 3456;
