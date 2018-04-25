param l {1..5};
param u {1..5};

var x {j in 1..5} >= l[j], <= u[j];

minimize obj: 
  exp( prod {j in 1..5} x[j] ) - 0.5*(x[1]^3 + x[2]^3 + 1)^2
  ;

subject to constr1: sum {j in 1..5} x[j]^2 = 10;
subject to constr2: x[2]*x[3] - 5*x[4]*x[5] = 0;
subject to constr3: x[1]^3 + x[2]^3 = -1;

data;

param l :=
  1  -2.3
  2  -2.3
  3  -3.2
  4  -3.2
  5  -3.2
  ;

param u :=
  1  2.3
  2  2.3
  3  3.2
  4  3.2
  5  3.2
  ;

let x[1] := -2;
let x[2] :=  2;
let x[3] :=  2;
let x[4] := -1;
let x[5] := -1;

#printf "optimal solution as starting point \n";
#let x[1] := -1.717142;
#let x[2] :=  1.159571;
#let x[3] :=  1.827248;
#let x[4] := -0.7636474;
#let x[5] := -0.7636390;

data;

solve;

display x;

display obj;

display obj - 0.0539498478;
