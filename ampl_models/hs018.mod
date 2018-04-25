var x {1..2};

minimize obj: 
  x[1]^2/100 + x[2]^2
  ;

subject to constr1: x[1]*x[2] >= 25;
subject to constr2: x[1]^2 + x[2]^2 >= 25;
subject to constr3: 2 <= x[1] <= 50;
subject to constr4: 0 <= x[2] <= 50;

let x[1] := 2;
let x[2] := 2;

#printf "optimal solution as starting point \n";
#let x[1] := 15.8114;
#let x[2] :=  1.58114;

solve;

display x;

display obj;

display obj - 5;
