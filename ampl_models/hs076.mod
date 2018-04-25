var x {j in 1..4} >= 0;

minimize obj: 
  x[1]^2 + 0.5*x[2]^2 + x[3]^2 + 0.5*x[4]^2 - x[1]*x[3] + x[3]*x[4]
  - x[1] - 3*x[2] + x[3] - x[4]
  ;

subject to constr1: x[1] + 2*x[2] + x[3] + x[4] <= 5;
subject to constr2: 3*x[1] + x[2] + 2*x[3] - x[4] <= 4;
subject to constr3: x[2] + 4*x[3] >= 1.5;

data;

let x[1] := 0.5;
let x[2] := 0.5;
let x[3] := 0.5;
let x[4] := 0.5;

#printf "optimal solution as starting point \n";
#let x[1] := 0.2727273;
#let x[2] := 2.090909;
#let x[3] :=  -0.26e-10;
#let x[4] := 0.5454545;

data;

solve;

display x;

display obj;

display obj + 4.681818181;
