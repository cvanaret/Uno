var x {1..4} >= 0;

minimize obj: 
  24.55*x[1] + 26.75*x[2] + 39*x[3] + 40.50*x[4]
  ;

subject to constr1: 2.3*x[1] + 5.6*x[2] + 11.1*x[3] + 1.3*x[4] >= 5;
subject to constr2: 12*x[1] + 11.9*x[2] + 41.8*x[3] + 52.1*x[4] 
	>= 21 + 
	   1.645*sqrt(0.28*x[1]^2 + 0.19*x[2]^2 + 20.5*x[3]^2 + 0.62*x[4]^2);
subject to constr3: sum {j in 1..4} x[j] = 1;

let x[1] := 1;
let x[2] := 1;
let x[3] := 1;
let x[4] := 1;

#printf "optimal solution as starting point \n";
#let x[1] := 0.6355216;
#let x[2] := -0.12e-11;
#let x[3] := 0.3127019;
#let x[4] := 0.5177655;

data;

solve;

display x;

display obj;

display obj - 29.894378;
