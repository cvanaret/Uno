param u {1..6}, default Infinity;
var x {j in 1..6} >= 0, <= u[j];

minimize obj: 
  x[1] + 2*x[2] + 4*x[5] + exp(x[1]*x[4]);
  ;

subject to constr1: x[1] + 2*x[2] + 5*x[5] = 6;
subject to constr2: x[1] + x[2] + x[3] = 3;
subject to constr3: x[4] + x[5] + x[6] = 2;
subject to constr4: x[1] + x[4] = 1;
subject to constr5: x[2] + x[5] = 2;
subject to constr6: x[3] + x[6] = 2;
#subject to constr7: x[1] <= 1;
#subject to constr8: x[4] <= 1;

data;

param u :=
    1  1
    4  1
    ;

let x[1] := 1;
let x[2] := 2;
let x[3] := 0;
let x[4] := 0;
let x[5] := 0;
let x[6] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 0;
#let x[2] := 4/3;
#let x[3] := 5/3;
#let x[4] := 1;
#let x[5] := 2/3;
#let x[6] := 1/3;

display obj;

solve;

display x;

display obj;

display obj - 19/3;
