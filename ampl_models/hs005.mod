var x {1..2};

minimize obj: 
  sin(x[1]+x[2]) + (x[1]-x[2])^2 - 1.5*x[1] + 2.5*x[2] + 1
  ;

subject to constr1: -1.5 <= x[1] <= 4;
subject to constr2: -3   <= x[2] <= 3;

let x[1] := 0;
let x[2] := 0;

solve;

display x;

display obj;

display obj + sqrt(3)/2 + 3.14159/3;
