var x {1..4} >= 0;

minimize obj: 
  x[1] - x[2] - x[3] - x[1]*x[3] + x[1]*x[4] + x[2]*x[3] - x[2]*x[4]
  ;

subject to constr1: x[1] + 2*x[2] <= 8;
subject to constr2: 4*x[1] + x[2] <= 12;
subject to constr3: 3*x[1] + 4*x[2] <= 12;
subject to constr4: 2*x[3] + x[4] <= 8;
subject to constr5: x[3] + 2*x[4] <= 8;
subject to constr6: x[3] + x[4] <= 5;

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
