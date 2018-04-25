var x {1..3} >= 1.0e-5;

minimize obj: 
  5*x[1] + 50000/x[1] + 20*x[2] + 72000/x[2] + 10*x[3] + 144000/x[3]
  ;

subject to constr1: 4/x[1] + 32/x[2] + 120/x[3] <= 1;

let x[1] := 1;
let x[2] := 1;
let x[3] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 108.7347175;
#let x[2] := 85.12613942;
#let x[3] := 204.3247078;

display obj;

solve;

display x;

display obj;

display obj - 6299.842428;
