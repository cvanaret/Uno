var x {1..5};

minimize obj: 
  (x[1]-1)^2 + (x[2]-x[3])^2 + (x[4]-x[5])^2
  ;

subject to constr1: sum {i in 1..5} x[i] = 5;
subject to constr2: x[3] - 2*(x[4]+x[5]) = -3;

let x[1] := 3;
let x[2] := 5;
let x[3] := -3;
let x[4] := 2;
let x[5] := -2;

display obj;

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj - 0;
