var x {1..3};


minimize obj: 
  4*x[1]^2 + 2*x[2]^2 + 2*x[3]^2 - 33*x[1] + 16*x[2] - 24*x[3]
  ;

subject to constr1: 3*x[1] - 2*x[2]^2 = 7;
subject to constr2: 4*x[1] -   x[3]^2 = 11;

let x[1] := 0;
let x[2] := 0;
let x[3] := 0;

#printf "optimal solution as starting point \n";
#let x[1] := 5.326770157;
#let x[2] := -2.118998639;
#let x[3] := 3.210464239;

display obj;

solve;

display x;

display obj;

display obj + 143.6461422;
