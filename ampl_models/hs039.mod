var x {1..4};

minimize obj: 
  -x[1]
  ;

subject to constr1: x[2] - x[1]^3 - x[3]^2 = 0;
subject to constr2: x[1]^2 - x[2] - x[4]^2 = 0;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
let x[4] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 1;
#let x[2] := 1;
#let x[3] := 0;
#let x[4] := 0;

display obj;

solve;

display x;

display obj;

display obj + 1;
