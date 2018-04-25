var x {1..5};

minimize obj: 
  (x[1]-x[2])^2 + (x[3]-1)^2 + (x[4]-1)^4 + (x[5]-1)^6
  ;

subject to constr1: x[1]^2*x[4] + sin(x[4] - x[5]) = 1;
subject to constr2: x[2] + x[3]^4*x[4]^2 = 2;

let x[1] := sqrt(2)/2;
let x[2] := 1.75;
let x[3] := 0.5;
let x[4] := 2;
let x[5] := 2;

display obj;

solve;

display x;

display obj;

display obj - 0;
