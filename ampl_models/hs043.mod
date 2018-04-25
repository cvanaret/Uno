var x {1..4};

minimize obj: 
  x[1]^2 + x[2]^2 + 2*x[3]^2 + x[4]^2 - 5*x[1] - 5*x[2] - 21*x[3] + 7*x[4]
  ;

subject to constr1: x[1]^2 +x[2]^2 +x[3]^2 +x[4]^2 +x[1] -x[2] +x[3] -x[4] <= 8;
subject to constr2: x[1]^2 +2*x[2]^2 +x[3]^2 +2*x[4]^2 -x[1] -x[4] <= 10;
subject to constr3: 2*x[1]^2 +x[2]^2 +x[3]^2 +2*x[1] -x[2] -x[4] <= 5;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;
let x[4] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 1;
#let x[3] := 2;
#let x[4] := -1;

display obj;

solve;

display x;

display obj;

display obj + 44;
