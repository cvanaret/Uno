param pi := 3.14159;

var x {1..2};

minimize obj: 
  sin(pi * x[1] / 12) * cos(pi * x[2] / 16);
  ;

subject to constr1: 4*x[1] - 3*x[2] = 0;

let x[1] := 0;
let x[2] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := -3;
#let x[2] := -4;

solve;

display x;

display obj;

display obj + 0.5;
