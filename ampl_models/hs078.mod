var x {1..5};

minimize obj: 
  prod {j in 1..5} x[j]
  ;

subject to constr1: sum {j in 1..5} x[j]^2 = 10;
subject to constr2: x[2]*x[3] - 5*x[4]*x[5] = 0;
subject to constr3: x[1]^3 + x[2]^3 = -1;

data;

let x[1] := -2;
let x[2] := 1.5;
let x[3] := 2;
let x[4] := -1;
let x[5] := -1;

#printf "optimal solution as starting point \n";
#let x[1] := -1.717142;
#let x[2] :=  1.595708;
#let x[3] :=  1.827248;
#let x[4] := -0.7636429;
#let x[5] := -0.7636435;

data;

solve;

display x;

display obj;

display obj + 2.91970041;
