var x {1..5};

minimize obj: 
  (x[1]-x[2])^2 + (x[2]-x[3])^2 + (x[3]-x[4])^4 + (x[4]-x[5])^2
  ;

subject to constr1: x[1] + 2*x[2] + 3*x[3] = 6;
subject to constr2: x[2] + 2*x[3] + 3*x[4] = 6;
subject to constr3: x[3] + 2*x[4] + 3*x[5] = 6;

let x[1] :=  35;
let x[2] := -31;
let x[3] :=  11;
let x[4] :=   5;
let x[5] :=  -5;

display obj;

solve;

display x;

display obj;

display obj - 0;
