var x {1..2};

minimize obj: 
  (x[1] - 5)^2 + x[2]^2 -25
  ;

subject to constr1: x[1]^2 <= x[2];

let x[1] := 4.9;
let x[2] := 0.1;

#printf "optimal solution as starting point \n";
#let x[1] := 1.23477;
#let x[2] := 1.52466;

solve;

display x;

display obj;

display obj + 8.498464223;
