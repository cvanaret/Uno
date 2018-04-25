var x {1..3} >= 0;

minimize obj: 
  -x[1]
  ;

subject to constr1: x[2] >= exp(x[1]);
subject to constr2: x[3] >= exp(x[2]);
subject to constr3: x[1] <= 100;
subject to constr4: x[2] <= 100;
subject to constr5: x[3] <= 10;

let x[1] := 0;
let x[2] := 1.05;
let x[3] := 2.9;

#printf "optimal solution as starting point \n";
#let x[1] := 0.83403;
#let x[2] := 2.30258;
#let x[3] := 10;

display obj;

solve;

display x;

display obj;

display obj + log(log(10));
