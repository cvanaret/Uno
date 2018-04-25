var x {1..2};

minimize obj: 
  log(1+x[1]^2) - x[2]
  ;

subject to constr: (1+x[1]^2)^2 + x[2]^2 = 4;

let x[1] := 2;
let x[2] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 1.73205;

solve;

display x;

display obj;

display obj + sqrt(3);
