var x {1..10} >= 2.001, <= 9.999;

minimize obj: 
  sum {j in 1..10} ( log(x[j]-2)^2 + log(10-x[j])^2 ) 
  - 
  (prod {j in 1..10} x[j])^0.2
  ;

let {j in 1..10} x[j] := 9;

solve;

display x;

display obj;

display obj + 45.77846971;
