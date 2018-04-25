var x {1..5};

minimize obj: 
  (x[1]-1)^2 + (x[1]-x[2])^2 + (x[2]-x[3])^2 + (x[3]-x[4])^4 + (x[4]-x[5])^4
  ;

subject to constr1: x[1] + x[2]^2 + x[3]^3 = 2 + 3*sqrt(2);
subject to constr2: x[2] - x[3]^2 + x[4]  = -2 + 2*sqrt(2);
subject to constr3: x[1]*x[5] = 2;

data;

let x[1] := 2;
let x[2] := 2;
let x[3] := 2;
let x[4] := 2;
let x[5] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 1.191127;
#let x[2] := 1.362603;
#let x[3] := 1.472818;
#let x[4] := 1.635017;
#let x[5] := 1.679081;

data;

solve;

display x;

display obj;

display obj - 0.0787768209;
