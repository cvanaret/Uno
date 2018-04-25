var x {1..4};

minimize obj: 
  -x[1]*x[2]*x[3]*x[4]
  ;

subject to constr1: x[1]^3 + x[2]^2 = 1;
subject to constr2: x[1]^2*x[4] - x[3] = 0;
subject to constr3: x[4]^2 - x[2] = 0;

let x[1] := 0.8;
let x[2] := 0.8;
let x[3] := 0.8;
let x[4] := 0.8;

#printf "optimal solution as starting point \n";
#let x[1] := 0.793701;
#let x[2] := 0.707107;
#let x[3] := 0.529732;
#let x[4] := 0.840896;

display obj;

solve;

display x;

display obj;

display obj + 0.25;
