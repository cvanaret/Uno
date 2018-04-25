var x {1..5};

minimize obj: 
  (x[1]-x[2])^2 + (x[2]-x[3])^3 + (x[3]-x[4])^4 + (x[4]-x[5])^4
  ;

subject to constr1: x[1] + x[2]^2 + x[3]^3 = 3;
subject to constr2: x[2] - x[3]^2 + x[4] = 1;
subject to constr3: x[1]*x[5] = 1;

let x[1] := 2;
let x[2] := sqrt(2);
let x[3] := -1;
let x[4] := 2-sqrt(2);
let x[5] := 1/2;

display obj;

solve;

display x;

display obj;

display obj - 0;
