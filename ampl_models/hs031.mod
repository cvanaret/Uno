var x {1..3};

minimize obj: 
  9*x[1]^2 + x[2]^2 + 9*x[3]^2
  ;

subject to constr1: x[1] * x[2] >= 1;
subject to constr2: -10 <= x[1] <= 10;
subject to constr3:   1 <= x[2] <= 10;
subject to constr4: -10 <= x[3] <= 1;

let x[1] :=  1;
let x[2] :=  1;
let x[3] :=  1;

#printf "optimal solution as starting point \n";
#let x[1] :=  0.57735;
#let x[2] :=  1.73205;
#let x[3] :=  0;

display obj;

solve;

display x;

display obj;

display obj - 6;
