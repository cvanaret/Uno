var x {1..7} >= 0;

param a := asin(sqrt(1/4.2));
param b := asin(sqrt(5/7.2));
param c := asin(sqrt(4/7));
param d := asin(sqrt(2/7));

minimize obj: 
  -x[1]*x[2]*x[3]
  ;

subject to constr1: x[1] - 4.2*sin(x[4])^2 = 0;
subject to constr2: x[2] - 4.2*sin(x[5])^2 = 0;
subject to constr3: x[3] - 4.2*sin(x[6])^2 = 0;
subject to constr4: x[1] + 2*x[2] + 2*x[3] - 7.2*sin(x[7])^2 = 0;

let x[1] := 1;
let x[2] := 1;
let x[3] := 1;
let x[4] := a;
let x[5] := a;
let x[6] := a;
let x[7] := b;

#printf "optimal solution as starting point \n";
#let x[1] := 2.4;
#let x[2] := 1.2;
#let x[3] := 1.2;
#let x[4] := c;
#let x[5] := d;
#let x[6] := d;
#let x[7] := 3.14159/2;

display obj;

solve;

display x;

display obj;

display obj + 3.456;
