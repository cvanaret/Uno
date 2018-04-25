var x {1..2};

minimize obj: 
  100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
  ;

subject to constr1: x[1]^2 + x[2] >= 0;
subject to constr2: x[1] + x[2]^2 >= 0;
subject to constr3: -1/2 <= x[1] <= 1/2;
subject to constr4: x[2] <= 1;

let x[1] := -2;
let x[2] :=  1;

#printf "optimal solution as starting point \n";
#let x[1] := 0.5;
#let x[2] := 0.25;

solve;

display x;

display obj;

display obj - 1/4;
