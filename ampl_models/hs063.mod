var x {1..3} >=0;

minimize obj: 
  1000 - x[1]^2 - 2*x[2]^2 - x[3]^2 - x[1]*x[2] - x[1]*x[3]
  ;

subject to constr1: 8*x[1] + 14*x[2] + 7*x[3] = 56;
subject to constr2: x[1]^2 +  x[2]^2 + x[3]^2 = 25;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 3.512118414;
#let x[2] := 0.2169881741;
#let x[3] := 3.552174034;

display obj;
display constr1.body, constr2.body;

solve;

display x;

display obj;

display obj - 961.7151721;
