var x {1..3};

minimize obj: 
  -x[1]*x[2]*x[3]
  ;

subject to constr1: x[1]^2 + 2*x[2]^2 + 4*x[3]^2 <= 48;

let x[1] :=  1;
let x[2] :=  1;
let x[3] :=  1;

#printf "optimal solution as starting point \n";
#let x[1] :=  4;
#let x[2] :=  2.82843;
#let x[3] :=  2;

display obj;

solve;

display x;

display obj;

display obj + 16*sqrt(2);
