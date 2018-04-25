var x {1..5} >= -10, <= 10;;

minimize obj: 
  (x[1]-x[2])^2 + (x[2]+x[3]-2)^2 + (x[4]-1)^2 + (x[5]-1)^2
  ;

subject to constr1: x[1] + 3*x[2] = 0;
subject to constr2: x[3] +   x[4] - 2*x[5] = 0;
subject to constr3: x[2] -   x[5] = 0;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
let x[4] := 2;
let x[5] := 2;

display obj;

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj - 176/43;
