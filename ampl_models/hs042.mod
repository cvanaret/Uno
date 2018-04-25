var x {1..4} >= 0;

minimize obj: 
  (x[1]-1)^2 + (x[2]-2)^2 + (x[3]-3)^2 + (x[4]-4)^2
  ;

subject to constr1: x[1] = 2;
subject to constr2: x[3]^2 + x[4]^2 = 2;

let x[1] := 1;
let x[2] := 1;
let x[3] := 1;
let x[4] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 2;
#let x[2] := 2;
#let x[3] := 0.848529;
#let x[4] := 1.13137;

display obj;

solve;

display x;

display obj;

display obj - 28 + 10*sqrt(2);
