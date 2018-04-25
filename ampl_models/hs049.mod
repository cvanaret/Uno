var x {1..5};

minimize obj: 
  (x[1]-x[2])^2 + (x[3]-1)^2 + (x[4]-1)^4 + (x[5]-1)^6
  ;

subject to constr1: sum {i in 1..4} x[i] + 3*x[4] = 7;
subject to constr2: x[3] + 5*x[5] = 6;

let x[1] := 10;
let x[2] :=  7;
let x[3] :=  2;
let x[4] := -3;
let x[5] := 0.8;

display obj;

solve;

display x;

display obj;

display obj - 0;
