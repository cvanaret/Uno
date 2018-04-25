var x {1..3};

minimize obj: 
  (x[1] + x[2])^2 + (x[2] + x[3])^2
  ;

subject to constr1: x[1] + 2*x[2] + 3*x[3] = 1;

let x[1] := -4;
let x[2] :=  1;
let x[3] :=  1;

#printf "optimal solution as starting point \n";
#let x[1] :=  0.5;
#let x[2] := -0.5;
#let x[3] :=  0.5;

display obj;

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj - 0;
