var x {1..2} <= 50, >= -50;

minimize obj: 
  x[1]^2 + x[2]^2
  ;

subject to constr1: x[1] + x[2] >= 1;
subject to constr2: x[1]^2 + x[2]^2 >= 1;
subject to constr3: 9*x[1]^2 + x[2]^2 >= 9;
subject to constr4: x[1]^2 - x[2] >= 0;
subject to constr5: x[2]^2 - x[1] >= 0;

let x[1] := 3;
let x[2] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 1;
#let x[2] := 1;

solve;

display x;

display obj;

display obj - 2;
