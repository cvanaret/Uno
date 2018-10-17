var x {1..2};
param l := 0;
param rho_ := 100;
var s >= 10;

minimize obj: 
	x[1]^2/100 + x[2]^2 - 100
   - l*(10*x[1] - x[2] - s)
   + rho_*(10*x[1] - x[2] - s)^2
  ;

subject to constr2: 2 <= x[1] <= 50;
subject to constr3: -50 <= x[2] <= 50;

let x[1] := -1;
let x[2] := -1;

#printf "optimal solution as starting point \n";
#let x[1] := 2.00265;
#let x[2] := 0;

solve;

display x;

display obj;

display obj + 99.96;
