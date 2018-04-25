var x {1..3};

param u {i in 1..99} := 25 + (-50*log(i/100))^(2/3);

minimize obj: 
  sum {i in 1..99} 
  (
      -i/100 + exp(-(u[i] - x[2])^x[3]/x[1])
  )^2
  ;

subject to constr1: 1/10 <= x[1] <= 100;
subject to constr2: 0 <= x[2] <= 25.6;
subject to constr3: 0 <= x[3] <= 5;

let x[1] := 100;
let x[2] := 12.5;
let x[3] :=  3;

#printf "optimal solution as starting point \n";
#let x[1] := ?
#let x[2] := ?
#let x[3] := ?

display obj;

solve;

display x;

display obj;

display obj - 0;
