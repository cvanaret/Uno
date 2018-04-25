var x {1..3} >= 0;

minimize obj: 
  (x[1] - 1)*(x[1] - 2)*(x[1] - 3) + x[3]
  ;

subject to constr1: x[1]^2 + x[2]^2 <= x[3]^2;
subject to constr2: x[1]^2 + x[2]^2 + x[3]^2 >= 4;
subject to constr3: x[3] <= 5;

let x[1] := 0;
let x[2] := 0;
let x[3] := 3;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := sqrt(2);
#let x[3] := sqrt(2);

display obj;

solve;

display x;

display obj;

display obj - sqrt(2) + 6;
