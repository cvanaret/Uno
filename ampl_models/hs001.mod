var x {1..2};

minimize obj: 
  100*(x[2] - x[1]^2)^2 + (1-x[1])^2;
  ;

subject to constr: -1.5 <= x[2];

let x[1] := -2;
let x[2] :=  1;

solve;

display x;

display obj;

display obj + 0.0;
