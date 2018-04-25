var x {1..5};

minimize obj: 
  (x[1]-1)^2 + (x[1] - x[2])^2 + (x[3]-1)^2 + (x[4]-1)^4 + (x[5]-1)^6
  ;

subject to constr1: x[1]^2*x[4] + sin(x[4]-x[5]) = 2*sqrt(2);
subject to constr2: x[2] + x[3]^4*x[4]^2 = 8 + sqrt(2);

data;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
let x[4] := 2;
let x[5] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 1.166172;
#let x[2] := 1.182111;
#let x[3] := 1.380257;
#let x[4] := 1.506036;
#let x[5] := 0.6109203;

data;

solve;

display x;

display obj;

display obj - 0.24150513;
