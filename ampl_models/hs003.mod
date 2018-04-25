var x {1..2};

minimize obj: 
  x[2] + 0.00001*(x[2]-x[1])^2
  ;

subject to constr: 0 <= x[2];

let x[1] := 10;
let x[2] :=  1;

option loqo_options $loqo_options" convex";

solve;

display x;

display obj;

display obj + 0.0;
