var x {1..2};

minimize obj: 
  (1-x[1])^2
  ;

subject to constr: 10*(x[2] - x[1]^2) = 0;

let x[1] := -1.2;
let x[2] :=  1;

#printf "optimal x as starting point \n";
#let x[1] := 1;
#let x[2] := 1;

solve;

display x;

display obj;

display obj + 0;
