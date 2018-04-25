var x {1..2};

minimize obj: 
  x[1] - x[2]
  ;

subject to constr1: -3*x[1]^2 + 2*x[1]*x[2] - x[2]^2 >= -1;

let x[1] := -10;
let x[2] :=  10;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 1;

#option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj + 1;
