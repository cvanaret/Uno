var x {1..3};

minimize obj: 
  (x[1] - x[2])^2 + (x[2] - x[3])^4
  ;

subject to constr1: (1 + x[2]^2)*x[1] + x[3]^4 = 3;

let x[1] := -2.6;
let x[2] :=  2;
let x[3] :=  2;

#printf "optimal solution as starting point \n";
#let x[1] :=  1;
#let x[2] :=  1;
#let x[3] :=  1;

display obj;

solve;

display x;

display obj;

display obj - 0;
