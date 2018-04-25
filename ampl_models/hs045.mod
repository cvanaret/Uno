var x {i in 1..5} >= 0, <= i;

minimize obj: 
  2 - (prod{i in 1..5} x[i])/120
  ;

#subject to constr1: x[1] + 2*x[2] <= 8;

let {i in 1..5} x[i] := 0;

display obj;

solve;

display x;

display obj;

display obj - 1;
