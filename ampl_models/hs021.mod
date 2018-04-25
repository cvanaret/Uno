var x {1..2};

minimize obj: 
  x[1]^2/100 + x[2]^2 - 100
  ;

subject to constr1: 10*x[1] - x[2] >= 10;
subject to constr2: 2 <= x[1] <= 50;
subject to constr3: -50 <= x[2] <= 50;

let x[1] := -1;
let x[2] := -1;

#printf "optimal solution as starting point \n";
#let x[1] := 2.00265;
#let x[2] := 0;

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj + 99.96;
