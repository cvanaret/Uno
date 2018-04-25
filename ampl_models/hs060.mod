var x {1..3} >= -10, <= 10;


minimize obj: 
  (x[1] - 1)^2 + (x[1] - x[2])^2 + (x[2] - x[3])^4
  ;

subject to constr1: x[1]*(1 + x[2]^2) + x[3]^4 = 4 + 3*sqrt(2);

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 1.104859024;
#let x[2] := 1.196674194;
#let x[3] := 1.535262257;

display obj;

solve;

display x;

display obj;

display obj - 0.03256820025;
