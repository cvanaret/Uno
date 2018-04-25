param l {1..5};
param u {1..5};

param a {1..12};

var x {j in 1..5} >= l[j], <= u[j];

minimize obj: 
  5.3578547*x[3]^2 + 0.8356891*x[1]*x[5] + 37.293239*x[1] - 40792.141
  ;

subject to constr1: 
    0 <= a[1] + a[2]*x[2]*x[5] + a[3]*x[1]*x[4] - a[4]*x[3]*x[5] <= 92;
subject to constr2: 
    0 <= a[5] + a[6]*x[2]*x[5] + a[7]*x[1]*x[2] + a[8]*x[3]^2 - 90 <= 20;
subject to constr3: 
    0 <= a[9] + a[10]*x[3]*x[5] + a[11]*x[1]*x[3] + a[12]*x[3]*x[4] -20 <= 5;

data;

param a :=
  1  85.334407
  2   0.0056858
  3   0.0006262
  4   0.0022053
  5  80.51249
  6   0.0071317
  7   0.0029955
  8   0.0021813
  9   9.300961
 10   0.0047026
 11   0.0012547
 12   0.0019085
  ;

param l :=
  1  78
  2  33
  3  27
  4  27
  5  27
  ;

param u :=
  1  102
  2  45
  3  45
  4  45
  5  45
  ;

let x[1] := 78;
let x[2] := 33;
let x[3] := 27;
let x[4] := 27;
let x[5] := 27;

#printf "optimal solution as starting point \n";
#let x[1] := 78;
#let x[2] := 33;
#let x[3] := 29.99526;
#let x[4] := 45;
#let x[5] := 36.77581;

solve;

display x;

display obj;

display obj + 30665.53867;
