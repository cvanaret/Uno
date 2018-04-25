var x {1..2};

minimize obj: 
  (x[1] - 2)^2 + (x[2] - 1)^2
  ;

subject to constr1: x[1] + x[2] <= 2;
subject to constr2: -x[1]^2 + x[2] >= 0;

let x[1] := 2;
let x[2] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 1;
#let x[2] := 1;

solve;

display x;

display obj;

display obj - 1;
