param u {1..2};

var x {j in 1..2} >= 0, <= u[j];


minimize obj: 
  -75.196 + 3.8112*x[1] + 0.0020567*x[1]^3 - 1.0345e-5*x[1]^4
  + 6.8306*x[2] - 0.030234*x[1]*x[2] + 1.28134e-3*x[2]*x[1]^2
  + 2.266e-7*x[1]^4*x[2] - 0.25645*x[2]^2 + 0.0034604*x[2]^3 - 1.3514e-5*x[2]^4
  + 28.106/(x[2] + 1) + 5.2375e-6*x[1]^2*x[2]^2 + 6.3e-8*x[1]^3*x[2]^2
  - 7e-10*x[1]^3*x[2]^3 - 3.405e-4*x[1]*x[2]^2 + 1.6638e-6*x[1]*x[2]^3
  + 2.8673*exp(0.0005*x[1]*x[2]) - 3.5256e-5*x[1]^3*x[2]
# the last term appears in CUTE but not in H&S
  -0.12694*x[1]^2
  ;

subject to constr1: x[1]*x[2] >= 700;
subject to constr2: x[2] - x[1]^2/125 >= 0;
subject to constr3: (x[2] - 50)^2 - 5*(x[1] - 55) >= 0;

data;

param u :=
  1  75
  2  65
  ;

let x[1] := 90;
let x[2] := 10;

#printf "optimal solution as starting point \n";
#let x[1] := 13.55010424;
#let x[2] := 51.66018129;

display obj;

solve;

display x;

display obj;

display obj + 7.804226324;
