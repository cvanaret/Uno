var x {1..2};

minimize obj: 
  100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
  ;

subject to constr1: x[1]*x[2] >= 1;
subject to constr2: x[1] + x[2]^2 >= 0;
subject to constr3: x[1] <= 1/2;

let x[1] := -2;
let x[2] :=  1;
