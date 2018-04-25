var x {1..2};

minimize obj: 
  (x[1] - 10)^3 + (x[2] - 20)^3
#  (x[1] - 10)* (x[1] - 10)* (x[1] - 10)
#  + 
#  (x[2] - 20)* (x[2] - 20)* (x[2] - 20)
  ;

subject to constr1: (x[1] - 5)^2 + (x[2] - 5)^2 >= 100;
subject to constr2: (x[2] - 5)^2 + (x[1] - 6)^2 <= 82.81;
subject to constr3: 13 <= x[1] <= 100;
subject to constr4: 0 <= x[2] <= 100;

let x[1] := 20.1; 
let x[2] :=  5.84;

#printf "optimal solution as starting point \n";
#let x[1] := 14.095;
#let x[2] :=  0.84296079;

solve;

display x;

display obj;

display obj + 6961.81381;
