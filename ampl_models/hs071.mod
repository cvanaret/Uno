var x {1..4} >= 1, <= 5;

minimize obj: 
  x[1]*x[4]*(x[1] + x[2] + x[3]) + x[3]
  ;

subject to constr1: prod {i in 1..4} x[i] >= 25;
subject to constr2: sum {i in 1..4} x[i]^2 = 40;

let x[1] := 1;
let x[2] := 5;
let x[3] := 5;
let x[4] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 1;
#let x[2] := 4.742994;
#let x[3] := 3.8211503;
#let x[4] := 1.3794082;

display obj;

solve;

display x;

display obj;

display obj - 17.0140173;
