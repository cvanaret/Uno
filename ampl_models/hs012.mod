var x {1..2};

minimize obj: 
  x[1]^2/2 + x[2]^2 - x[1]*x[2] - 7*x[1] - 7*x[2]
  ;

subject to constr1: 4*x[1]^2 + x[2]^2 <= 25;

let x[1] := 0;
let x[2] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 2;
#let x[2] := 3;

solve;

display x;

display obj;

display obj + 30;
