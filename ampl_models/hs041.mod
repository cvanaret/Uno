var x {1..4} >= 0;

minimize obj: 
  2-x[1]*x[2]*x[3]
  ;

subject to constr1: x[1] + 2*x[2] + 2*x[3] - x[4] = 0;
subject to constr2: x[1] <= 1;
subject to constr3: x[2] <= 1;
subject to constr4: x[3] <= 1;
subject to constr5: x[4] <= 2;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
let x[4] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 2/3;
#let x[2] := 1/3;
#let x[3] := 1/3;
#let x[4] := 2;

display obj;

solve;

display x;

display obj;

display obj - 52/27;
